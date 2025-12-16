#!/usr/bin/env python3
"""
Automated generation of barcode fasta files.

For each barcode type with 'generate' option set, generates a fasta file containing
the most frequently appearing barcode sequences from the provided BAM file.
Sequences are chosen based on frequency, with only the most frequent sequences
written to the file. If hamming_distance option is set, all sequences written
will be of hamming distance > specified distance from each other. If the barcode
fasta file already exists, a new one will not be generated.

The script can run standalone or as part of a Snakemake pipeline.
"""

import os
import sys
import argparse
import subprocess
import tempfile
from pathlib import Path
from collections import Counter

import pandas as pd
import pysam
from Bio import Seq, SeqIO

# Import from common and demux_refactored
from common import load_csv_as_dict, validate_bc
from demux import build_reference_alignment, find_barcode_position_in_alignment, hamming_distance


def find_barcodes_in_bam(bam_path, reference_fasta, barcode_info):
    """
    Find all barcodes of specified types in BAM file.

    Args:
        bam_path: Path to input BAM file
        reference_fasta: Path to reference FASTA file
        barcode_info: Dictionary of barcode information

    Returns:
        tuple: (barcode_counts, barcode_qualities, has_quality_scores)
            - barcode_counts: {barcode_type: {sequence: count}}
            - barcode_qualities: {barcode_type: {sequence: quality_string}}
            - has_quality_scores: bool indicating if BAM has quality scores
    """
    # Load all reference sequences
    references = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(reference_fasta, 'fasta')}

    # Initialize barcode counts and qualities for types that need generation
    barcode_counts = {}
    barcode_qualities = {}
    for barcode_type, info in barcode_info.items():
        if info.get('generate', False):
            barcode_counts[barcode_type] = Counter()
            barcode_qualities[barcode_type] = {}

    if not barcode_counts:
        return {}, {}, False

    # Open BAM file and process reads
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    has_quality_scores = None

    # Process reads from all references
    for ref_name, reference_seq in references.items():
        for bam_entry in bam_file.fetch(ref_name):
            # Check if quality scores exist (only need to check once)
            if has_quality_scores is None:
                has_quality_scores = bam_entry.query_qualities is not None

            ref_alignment = build_reference_alignment(bam_entry, reference_seq)

            for barcode_type, counts in barcode_counts.items():
                info = barcode_info[barcode_type]
                context = info['context'].upper()

                start, end = find_barcode_position_in_alignment(ref_alignment, context)

                if isinstance(start, int):
                    # Extract barcode from query sequence
                    barcode_seq = bam_entry.query_alignment_sequence[start:end]

                    # Apply reverse complement if needed
                    if info.get('reverse_complement', False):
                        barcode_seq = str(Seq.reverse_complement(barcode_seq))

                    # Skip if contains N
                    if 'N' not in barcode_seq:
                        counts[barcode_seq] += 1

                        # Store quality scores (keep first occurrence)
                        if barcode_seq not in barcode_qualities[barcode_type]:
                            if has_quality_scores:
                                qual_scores = bam_entry.query_qualities[start:end]
                                # Convert to ASCII phred+33 format
                                qual_string = ''.join(chr(q + 33) for q in qual_scores)
                            else:
                                # Fake quality scores (all Q40)
                                qual_string = 'I' * len(barcode_seq)
                            barcode_qualities[barcode_type][barcode_seq] = qual_string

    bam_file.close()

    return barcode_counts, barcode_qualities, has_quality_scores if has_quality_scores is not None else False


def filter_barcodes_by_hamming_distance(barcode_df, min_hamming_distance, quality_dict):
    """
    Filter barcodes to ensure minimum hamming distance between selected sequences.
    Uses UMICollapse for efficient clustering.

    Args:
        barcode_df: DataFrame with columns ['barcode', 'count'] sorted by count descending
        min_hamming_distance: Minimum hamming distance required
        quality_dict: Dict of {sequence: quality_string} (real or fake quality scores)

    Returns:
        list: Selected barcode sequences
    """
    if min_hamming_distance == 0 or len(barcode_df) == 0:
        return barcode_df['barcode'].tolist()

    # Create tmp directory
    os.makedirs('tmp', exist_ok=True)

    # Create temp files
    temp_input = tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False, dir='tmp')
    temp_output = tempfile.NamedTemporaryFile(mode='r', suffix='.fastq', delete=False, dir='tmp')
    temp_input_path = temp_input.name
    temp_output_path = temp_output.name
    temp_input.close()
    temp_output.close()

    # Write barcodes to FASTQ
    with open(temp_input_path, 'w') as f:
        for idx, row in barcode_df.iterrows():
            barcode = row['barcode']
            count = row['count']
            f.write(f'@bc{idx}_count{count}\n')
            f.write(f'{barcode}\n')
            f.write('+\n')
            f.write(f'{quality_dict[barcode]}\n')

    # Run UMICollapse
    cmd = [
        'umicollapse',
        'fastq',
        '-i', temp_input_path,
        '-o', temp_output_path,
        '-k', str(min_hamming_distance)
    ]
    subprocess.run(cmd, check=True, capture_output=True)

    # Parse output to get selected barcodes
    selected_barcodes = []
    for record in SeqIO.parse(temp_output_path, 'fastq'):
        selected_barcodes.append(str(record.seq))

    # Cleanup temp files
    os.unlink(temp_input_path)
    os.unlink(temp_output_path)

    # Remove tmp directory if empty
    try:
        os.rmdir('tmp')
    except OSError:
        pass  # Directory not empty or doesn't exist

    return selected_barcodes


def write_barcode_fasta(barcode_counts, barcode_info, barcode_qualities):
    """
    Write barcode FASTA files for each barcode type.

    Args:
        barcode_counts: Dictionary of barcode counts {type: {seq: count}}
        barcode_info: Dictionary of barcode information
        barcode_qualities: Dictionary of barcode quality scores {type: {seq: quality_string}}

    Returns:
        None
    """
    if not barcode_counts:
        return

    for barcode_type, counts in barcode_counts.items():
        if not counts:
            continue

        info = barcode_info[barcode_type]
        fasta_path = info['fasta']
        max_barcodes = info.get('generate', 'all')
        hamming_dist = info.get('hamming_distance', 0)

        # Convert counts to sorted DataFrame
        barcode_list = [(seq, count) for seq, count in counts.items()]
        df = pd.DataFrame(barcode_list, columns=['barcode', 'count'])
        df = df.sort_values('count', ascending=False).reset_index(drop=True)

        # Filter by hamming distance if specified
        if hamming_dist > 0:
            quality_dict = barcode_qualities.get(barcode_type, {})
            selected_barcodes = filter_barcodes_by_hamming_distance(df, hamming_dist, quality_dict)
            # Filter DataFrame to only selected barcodes
            df = df[df['barcode'].isin(selected_barcodes)]

        # Create output directory if needed
        output_dir = os.path.dirname(fasta_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # Write FASTA file
        with open(fasta_path, 'w') as f:
            count = 0
            for _, row in df.iterrows():
                count += 1
                f.write(f'>bc{count}\n')
                f.write(f'{row.barcode}\n')

                # Stop if reached max barcodes
                if max_barcodes != 'all' and count >= max_barcodes:
                    break


def main():
    # Parse arguments
    if 'snakemake' in globals():
        # Running as snakemake script
        # bam and references are lists (one per tag)
        bam_files = snakemake.input.bam if isinstance(snakemake.input.bam, list) else [snakemake.input.bam]
        references = snakemake.params.references if isinstance(snakemake.params.references, list) else [snakemake.params.references]
        tags = snakemake.params.tags if isinstance(snakemake.params.tags, list) else [snakemake.params.tags]

        args = type('Args', (), {
            'bam': bam_files,
            'reference': references,
            'barcode_info': snakemake.input.barcode_info,
            'tags': tags,
            'output_fasta': snakemake.output.fasta
        })()
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--bam', required=True, nargs='+', help='Input BAM file(s)')
        parser.add_argument('--reference', required=True, nargs='+', help='Reference FASTA file(s) (one per BAM)')
        parser.add_argument('--barcode_info', required=True, help='Barcode info CSV file')
        parser.add_argument('--output_fasta', required=True, help='Output barcode FASTA file')
        args = parser.parse_args()
        args.tags = [f'tag{i}' for i in range(len(args.bam))]  # Dummy tags for standalone

    if args.barcode_info:
        # Load barcode information (use first tag for validation, but doesn't really matter)
        barcode_info = load_csv_as_dict(
            args.barcode_info,
            required=['barcode_name', 'context', 'fasta'],
            tag=args.tags[0],
            defaults={'reverse_complement': False, 'generate': False, 'hamming_distance': 0}
        )

        # Validate barcode info
        metadata_dir = Path(args.barcode_info).parts[0]
        for bc in barcode_info:
            barcode_info, _ = validate_bc(barcode_info, bc, metadata_dir=metadata_dir)

        # Aggregate barcode counts and qualities across all BAM files
        all_barcode_counts = {}
        all_barcode_qualities = {}
        for bam_file, reference in zip(args.bam, args.reference):
            barcode_counts, barcode_qualities, has_quality = find_barcodes_in_bam(bam_file, reference, barcode_info)

            # Merge counts from this BAM into all_barcode_counts
            for barcode_type, type_counts in barcode_counts.items():
                if barcode_type not in all_barcode_counts:
                    all_barcode_counts[barcode_type] = Counter()
                    all_barcode_qualities[barcode_type] = {}
                all_barcode_counts[barcode_type].update(type_counts)

                # Merge qualities (keep first occurrence)
                for seq, qual in barcode_qualities[barcode_type].items():
                    if seq not in all_barcode_qualities[barcode_type]:
                        all_barcode_qualities[barcode_type][seq] = qual

        # Write barcode FASTA files
        write_barcode_fasta(all_barcode_counts, barcode_info, all_barcode_qualities)



if __name__ == '__main__':
    main()