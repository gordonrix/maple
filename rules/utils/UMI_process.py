#!/usr/bin/env python3
"""
UMI_process.py

UMI processing script for the maple pipeline.

Given an aligned BAM file containing reads aligned to multiple reference sequences:
1. Extracts UMI sequences from each read based on UMI contexts in the reference it aligned to
2. Groups reads by (reference_name, UMI_sequence) with mismatch tolerance using UMICollapse for graph-based clustering
3. Filters groups by minimum/maximum read count
4. Splits groups into batches for consensus generation
5. Outputs batch FASTA files and log files

Supports both snakemake mode and standalone command-line mode.

Author: Gordon Rix
"""

import argparse
import pandas as pd
import numpy as np
import pysam
import os
import sys
import shutil
import subprocess
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import bisect

# Import from existing maple utilities
from demux import build_reference_alignment, find_barcode_position_in_alignment

# Check if running in snakemake mode
try:
    snakemake
    snakemake_mode = True
except NameError:
    snakemake_mode = False

# Set up arguments for snakemake or command line
if snakemake_mode:
    class Args: pass
    args = Args()
    args.input_bam = snakemake.input.bam
    args.references = snakemake.input.references
    args.output_bams = snakemake.output.bams
    args.extract_log = snakemake.output.extract_log
    args.groups_log = snakemake.output.groups_log
    args.UMI_contexts = snakemake.params.UMI_contexts
    args.UMI_mismatches = snakemake.params.UMI_mismatches
    args.batches = snakemake.params.batches
    args.minimum = snakemake.params.minimum
    args.maximum = snakemake.params.maximum
else:
    parser = argparse.ArgumentParser(description="Unified UMI processing: extract, group, filter, and batch")
    parser.add_argument("--input_bam", required=True, help="Input aligned BAM file")
    parser.add_argument("--references", required=True, help="FASTA file containing reference sequences")
    parser.add_argument("--output_dir", required=True, help="Output directory for batch BAM files")
    parser.add_argument("--extract_log", required=True, help="Output CSV log for UMI extraction")
    parser.add_argument("--groups_log", required=True, help="Output CSV log for UMI groups")
    parser.add_argument("--UMI_contexts", required=True, help="Comma-separated UMI contexts")
    parser.add_argument("--UMI_mismatches", type=int, required=True, help="Maximum UMI mismatches for grouping")
    parser.add_argument("--batches", type=int, required=True, help="Number of batches to split into")
    parser.add_argument("--minimum", type=int, required=True, help="Minimum reads per UMI group")
    parser.add_argument("--maximum", type=int, required=True, help="Maximum reads per UMI group")
    args = parser.parse_args()

# Load references from file
references = {}
with pysam.FastxFile(args.references) as fh:
    for entry in fh:
        references[entry.name] = entry.sequence.upper()


def extract_umi(aligned_ref, bam_entry, UMI_contexts):
    """
    Extract concatenated UMI sequence from a read.

    Args:
        aligned_ref (str): Aligned reference sequence
        bam_entry: pysam AlignedSegment
        UMI_contexts (list): List of UMI context strings

    Returns:
        tuple: (umi_string, failure_flags) where umi_string is concatenated UMI or None,
               and failure_flags is a list of bools indicating failure per context
    """
    umi_string = ''
    failure_flags = []

    for context in UMI_contexts:
        result = find_barcode_position_in_alignment(aligned_ref, context)

        # Check if extraction failed
        if result[0] is None:
            failure_flags.append(True)
            return None, failure_flags

        start_pos, end_pos = result
        umi_string += bam_entry.query_alignment_sequence[start_pos:end_pos]
        failure_flags.append(False)

    return umi_string, failure_flags


def extract_umis_from_bam(bam_path, references, UMI_contexts):
    """
    Extract UMI sequences from aligned BAM file.

    For each read in the BAM:
    1. Identify which reference it aligned to (from RNAME field)
    2. Reconstruct the aligned reference sequence using CIGAR
    3. Extract UMI sequences based on UMI contexts in that reference
    4. Store read data grouped by (reference_name, UMI_sequence)

    Args:
        bam_path (str): Path to aligned BAM file
        references (dict): Dictionary of reference sequences {ref_name: {'alignment_seq': seq, ...}}
        UMI_contexts (list): List of UMI context strings (e.g., ['NNNYRNNN', 'NNNYRNNN'])

    Returns:
        tuple: (umi_groups, extract_log_entries)
            - umi_groups: dict {(ref_name, umi_seq): [ReadData, ...]}
            - extract_log_entries: list of [read_id, umi, ref_name, success, failure, umi_1_failure, ...]
    """
    umi_groups = defaultdict(list)
    extract_log_entries = []

    bam_in = pysam.AlignmentFile(bam_path, 'rb')

    for read in bam_in.fetch():
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        ref_name = read.reference_name

        # Skip if reference not in our references dict
        if ref_name not in references:
            continue

        reference_seq = references[ref_name]

        # Reconstruct aligned reference
        aligned_ref = build_reference_alignment(read, reference_seq)

        # Extract UMI
        umi, failure_flags = extract_umi(aligned_ref, read, UMI_contexts)

        if umi is None:
            # Log failure
            extract_log_entries.append([read.query_name, '', ref_name, 0, 1] + [int(f) for f in failure_flags])
            continue

        # Create ReadData object
        strand = '-' if read.is_reverse else '+'
        quality = read.query_qualities if read.query_qualities else [40] * len(read.query_sequence)
        read_data = ReadData(
            read_id=read.query_name,
            sequence=read.query_sequence,
            quality=quality,
            strand=strand,
            reference_name=ref_name,
            umi=umi
        )

        # Add to groups
        umi_groups[(ref_name, umi)].append(read_data)

        # Log success
        extract_log_entries.append([read.query_name, umi, ref_name, 1, 0] + [int(f) for f in failure_flags])

    bam_in.close()

    return umi_groups, extract_log_entries


def group_umis(umi_groups, UMI_mismatches, output_dir):
    """
    Group UMIs with similar sequences together (allowing mismatches) using UMICollapse.

    For each reference separately:
    1. Write UMI FASTQ to temp file
    2. Run UMICollapse to cluster UMIs with mismatch tolerance
    3. Parse UMICollapse output to get unique_id assignments
    4. Reorganize groups to use (ref_name, unique_id) as key
    5. Clean up temp directory

    Args:
        umi_groups (dict): {(ref_name, umi_seq): [ReadData, ...]}
        UMI_mismatches (int): Maximum mismatches allowed when grouping UMIs
        output_dir (str): Output directory (temp files go in output_dir/tmp{random}/)

    Returns:
        tuple: (grouped_umi_groups, groups_log_entries)
            - grouped_umi_groups: dict {(ref_name, unique_id): [ReadData, ...]}
            - groups_log_entries: list of [read_id, umi, ref_name, umi_count, final_umi_count, unique_id]
    """
    # Generate random directory name that doesn't exist
    while True:
        random_id = random.randint(100000000, 999999999)
        tmp_dir = os.path.join(output_dir, f'tmp{random_id}')
        if not os.path.exists(tmp_dir):
            break
    os.makedirs(tmp_dir)

    grouped_umi_groups = {}
    groups_log_entries = []

    # Get unique references
    references = set(ref_name for ref_name, _ in umi_groups.keys())

    for ref_name in references:
        # Get all UMI groups for this reference
        ref_umi_groups = {(rn, umi): reads for (rn, umi), reads in umi_groups.items() if rn == ref_name}

        if not ref_umi_groups:
            continue

        # Create temp FASTQ files
        tmp_umi_in = os.path.join(tmp_dir, f"{ref_name}_umis.fastq")
        tmp_umi_out = os.path.join(tmp_dir, f"{ref_name}_umis_grouped.fastq")

        # Write UMIs to FASTQ (one entry per read with UMI)
        with open(tmp_umi_in, 'w') as f:
            for (rn, umi), read_list in ref_umi_groups.items():
                for read in read_list:
                    # Write FASTQ record: read_id contains original read ID
                    f.write(f"@{read.read_id}\n")
                    f.write(f"{umi}\n")
                    f.write(f"+\n")
                    f.write(f"{'I' * len(umi)}\n")  # Dummy quality scores

        # Run UMICollapse
        print(f"Running UMI collapse on reference sequence {ref_name}")
        subprocess.run([
            "umicollapse", "fastq",
            "-i", tmp_umi_in,
            "-o", tmp_umi_out,
            "--tag",
            "-k", str(UMI_mismatches)
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Parse UMICollapse output to get cluster assignments
        umi_to_cluster = {}  # {(read_id, umi): cluster_id}

        for record in SeqIO.parse(tmp_umi_out, "fastq"):
            # Parse description for UMICollapse tags
            desc_dict = {}
            for item in record.description.split()[1:]:
                if '=' in item:
                    key, value = item.split('=', 1)
                    desc_dict[key] = value

            read_id = record.id
            umi = str(record.seq)
            cluster_id = int(desc_dict['cluster_id'])
            umi_count = desc_dict.get('same_umi', pd.NA)
            final_umi_count = desc_dict.get('cluster_size', pd.NA)

            # Store cluster assignment
            umi_to_cluster[(read_id, umi)] = cluster_id

            # Add to groups log
            groups_log_entries.append([read_id, umi, ref_name, umi_count, final_umi_count, cluster_id])

        # Reorganize reads by (ref_name, cluster_id)
        for (_, umi), read_list in ref_umi_groups.items():
            for read in read_list:
                cluster_id = umi_to_cluster.get((read.read_id, umi))
                if cluster_id is not None:
                    key = (ref_name, cluster_id)
                    if key not in grouped_umi_groups:
                        grouped_umi_groups[key] = []
                    grouped_umi_groups[key].append(read)

    # Clean up temp directory
    shutil.rmtree(tmp_dir)

    return grouped_umi_groups, groups_log_entries


def filter_and_batch_groups(umi_groups, bam_in, output_dir, minimum, maximum, batches):
    """
    Filter UMI groups by size and split into BAM batch files for consensus.

    Writes BAM files preserving alignment information from the input BAM.

    Args:
        umi_groups (dict): {(ref_name, unique_id): [ReadData, ...]}
        bam_in (pysam.AlignmentFile): Input BAM file (must be open)
        output_dir (str): Directory for output batch BAM files
        minimum (int): Minimum reads required per UMI group
        maximum (int): Maximum reads to keep per UMI group
        batches (int): Number of batch files to create
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Create dict: read_id → BAM record for fast lookup
    print("Building read ID to BAM record mapping...")
    read_to_record = {}
    for read in bam_in.fetch():
        if not read.is_unmapped:
            read_to_record[read.query_name] = read
    print(f"Loaded {len(read_to_record)} BAM records")

    # Open all batch BAM file handles
    batch_handles = {}
    for batch_num in range(batches):
        batch_path = os.path.join(output_dir, f'batch{batch_num}.bam')
        # Create BAM with same header as input
        batch_handles[batch_num] = pysam.AlignmentFile(batch_path, 'wb', template=bam_in)

    # Process each UMI group (same filtering logic as FASTA version)
    groups_written = 0
    reads_written = 0
    for (ref_name, unique_id), read_list in umi_groups.items():
        # Filter by minimum
        if len(read_list) < minimum:
            continue

        # Sort by quality (ascending, so lowest quality first)
        read_list.sort(key=lambda r: r.mean_quality)

        # Downsample to maximum if needed
        if len(read_list) > maximum:
            # Build parallel strand list
            group_strand_dict = [1 if r.strand == '+' else -1 for r in read_list]

            # Remove lowest quality reads while balancing strand bias
            while len(read_list) > maximum:
                strand_bias = sum(group_strand_dict)
                if len(set(group_strand_dict)) == 2:
                    # Both strands present - remove lowest quality of over-represented strand
                    if strand_bias < 0:
                        remove_idx = group_strand_dict.index(-1)
                    else:
                        remove_idx = group_strand_dict.index(1)
                else:
                    # Only one strand - remove lowest quality
                    remove_idx = 0

                read_list.pop(remove_idx)
                group_strand_dict.pop(remove_idx)

        # Assign to batch
        batch_num = unique_id % batches

        # Write BAM records to batch file
        for read in read_list:
            if read.read_id in read_to_record:
                bam_record = read_to_record[read.read_id]
                # Update read name to include UMI info
                bam_record.query_name = f'UMI-{unique_id}|ref-{ref_name}|{read.read_id}'
                batch_handles[batch_num].write(bam_record)
                reads_written += 1
            else:
                print(f"Warning: Read {read.read_id} not found in BAM file")
        groups_written += 1

    # Close all batch BAM files
    for handle in batch_handles.values():
        handle.close()

    # Sort and index all batch BAM files
    for batch_num in range(batches):
        batch_path = os.path.join(output_dir, f'batch{batch_num}.bam')
        sorted_path = os.path.join(output_dir, f'batch{batch_num}.sorted.bam')

        # Sort the BAM file
        pysam.sort("-o", sorted_path, batch_path)

        # Replace unsorted with sorted
        os.remove(batch_path)
        os.rename(sorted_path, batch_path)

        # Index the sorted BAM file
        pysam.index(batch_path)

    print(f"Wrote {reads_written} reads from {groups_written} UMI groups to {batches} batch BAM files")


def write_extract_log(extract_log_entries, log_path, num_umi_contexts):
    """
    Write UMI extraction log to CSV.

    Args:
        extract_log_entries (list): List of log entries
        log_path (str): Output CSV path
        num_umi_contexts (int): Number of UMI contexts (for column headers)
    """
    columns = ['read_id', 'umi', 'reference_name', 'success', 'failure'] + [f'umi_{i+1}_failure' for i in range(num_umi_contexts)]
    df = pd.DataFrame(extract_log_entries, columns=columns)
    df.to_csv(log_path, index=False)


def write_groups_log(groups_log_entries, log_path):
    """
    Write UMI groups log to CSV.

    Follows the same efficient pattern as UMI_groups_log.py but adds reference_name column.
    Expected entries: [read_id, umi, reference_name, umi_count, final_umi_count, unique_id]
    where umi_count and final_umi_count are only populated for consensus reads (NA for others).

    Args:
        groups_log_entries (list): List of log entries from UMICollapse parsing
        log_path (str): Output CSV path
    """
    columns = ['read_id', 'umi', 'reference_name', 'umi_count', 'final_umi_count', 'unique_id']
    df = pd.DataFrame(groups_log_entries, columns=columns)

    # Helper columns for sorting (to bring non-NA values to the front)
    df["is_na_fumi"] = df["final_umi_count"].isna().astype(int)
    df["is_na_umi"] = df["umi_count"].isna().astype(int)

    # Sort to enable consensus filling within groups
    df = df.sort_values(
        by=["reference_name", "unique_id", "is_na_fumi", "final_umi_count", "umi", "is_na_umi", "umi_count"],
        ascending=[True, True, True, False, True, True, False],
        na_position="last",
        kind="stable"
    ).reset_index(drop=True)
    df = df.drop(columns=["is_na_fumi", "is_na_umi"])

    # Forward fill consensus values within each (reference_name, unique_id) group
    df["final_umi_count"] = df.groupby(['reference_name', 'unique_id'])["final_umi_count"].ffill().bfill()
    df["umi_count"] = df.groupby(['reference_name', 'unique_id'])["umi_count"].ffill().bfill()

    # Convert to numeric for proper sorting
    df["final_umi_count"] = pd.to_numeric(df["final_umi_count"])
    df["umi_count"] = pd.to_numeric(df["umi_count"])
    df["unique_id"] = pd.to_numeric(df["unique_id"])

    # Final sort
    df = df.sort_values(by=['reference_name', 'final_umi_count', 'unique_id', 'umi_count'],
                        ascending=[True, False, True, False])
    df.to_csv(log_path, index=False)


class ReadData:
    """
    Container for read data needed for consensus generation.
    """
    def __init__(self, read_id, sequence, quality, strand, reference_name, umi):
        self.read_id = read_id
        self.sequence = sequence
        self.quality = quality  # List of quality scores
        self.strand = strand  # '+' or '-'
        self.reference_name = reference_name
        self.umi = umi
        self.mean_quality = np.mean(quality) if quality else 0


def main():
    """
    Main entry point for UMI processing.
    """
    # Get output directory from first bam path
    if snakemake_mode:
        output_dir = os.path.dirname(args.output_bams[0])
    else:
        output_dir = args.output_dir

    # Normalize UMI contexts
    UMI_contexts = [c.upper() for c in args.UMI_contexts]

    # Step 1: Extract UMIs from BAM file
    print("Step 1: Extracting UMIs from BAM file...")
    umi_groups, extract_log_entries = extract_umis_from_bam(
        args.input_bam,
        references,
        UMI_contexts
    )

    # Step 2: Group UMIs with mismatch tolerance using UMICollapse
    print("Step 2: Grouping UMIs with mismatch tolerance...")
    umi_groups, groups_log_entries = group_umis(
        umi_groups,
        args.UMI_mismatches,
        output_dir
    )

    # Step 3: Filter and batch groups
    print("Step 3: Filtering and batching groups...")
    # Open BAM file for reading to get records
    bam_in = pysam.AlignmentFile(args.input_bam, 'rb')
    filter_and_batch_groups(
        umi_groups,
        bam_in,
        output_dir,
        args.minimum,
        args.maximum,
        args.batches
    )
    bam_in.close()

    # Step 4: Write logs
    print("Step 4: Writing logs...")
    write_extract_log(extract_log_entries, args.extract_log, len(UMI_contexts))
    write_groups_log(groups_log_entries, args.groups_log)

    print("UMI processing complete!")


if __name__ == '__main__':
    main()
