#!/usr/bin/env python3
"""
Consensus Read Validation with Aligned Protein Track

This script is a tool that can be used with the Maple pipeline to retrieve and visualize
the original raw sequencing reads that were fed to medaka to generate a consensus sequence.
It processes multiple samples from a CSV input file and outputs one FASTA and PNG file for each row
"""

import pandas as pd
import os
import sys
import argparse
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from pathlib import Path


def reconstruct_alignment_from_bam(read, ref_len):
    """
    Reconstructs the read sequence anchored to the reference coordinates.

    Args:
        read: pysam AlignedSegment object
        ref_len: Length of the reference sequence

    Returns:
        tuple: (aligned_sequence_string, set_of_insertion_locations)
    """
    aligned_seq = [' '] * ref_len
    insertion_locs = []
    aligned_pairs = read.get_aligned_pairs(matches_only=False)
    read_seq = read.query_sequence

    last_ref_pos = -1
    for query_pos, ref_pos in aligned_pairs:
        if ref_pos is not None:
            if ref_pos < ref_len:
                aligned_seq[ref_pos] = read_seq[query_pos] if query_pos is not None else '-'
                last_ref_pos = ref_pos
        else:
            if last_ref_pos != -1 and last_ref_pos < ref_len:
                insertion_locs.append(last_ref_pos)

    return "".join(aligned_seq), set(insertion_locs)


def calculate_protein_mapping(ref_aln_str, ref_coding_str):
    """
    Aligns the coding sequence to the alignment reference to determine
    where to plot amino acids and their indices.

    This function handles cases where the coding sequence may be a substring
    of or differ from the alignment reference, using local alignment to find
    the correct reading frame.

    Args:
        ref_aln_str: The reference sequence used for alignment display
        ref_coding_str: The coding sequence (CDS) for protein translation

    Returns:
        list: List of dicts with keys 'plot_idx' (position in alignment),
              'aa' (amino acid character), and 'aa_idx' (1-based amino acid index)
    """
    # Translate the full coding sequence
    full_protein = Seq(ref_coding_str).translate()

    # Align coding DNA to reference DNA using local alignment
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -2

    alignments = aligner.align(ref_aln_str, ref_coding_str)
    best_aln = alignments[0]

    # Get coordinates of the match in both sequences
    target_start, target_end = best_aln.aligned[0][0]
    query_start, query_end = best_aln.aligned[1][0]

    mapping = []

    # Iterate through the aligned region and map codons to amino acids
    for plot_idx in range(target_start, target_end):
        coding_idx = query_start + (plot_idx - target_start)

        # Plot amino acid at the start of each codon
        if coding_idx % 3 == 0:
            aa_index = coding_idx // 3
            if aa_index < len(full_protein):
                aa_char = full_protein[aa_index]
                mapping.append({
                    'plot_idx': plot_idx,
                    'aa': aa_char,
                    'aa_idx': aa_index + 1  # 1-based index for display
                })

    return mapping


def get_original_sequences_aligned(tag, barcode, base_maple_dir, label_barcode=None, seq_id=None):
    """
    Extracts and aligns original reads for a specific consensus sequence.

    Args:
        tag: Sample tag identifier
        barcode: Barcode group identifier
        base_maple_dir: Base directory containing MAPLE pipeline outputs
        label_barcode: Specific barcode label to find (optional if seq_id provided)
        seq_id: Direct sequence ID (optional, takes precedence over label_barcode)

    Returns:
        tuple: (aligned_reads_data, reference_alignment_seq, reference_coding_seq, metadata_dict)

    Raises:
        FileNotFoundError: If required metadata files are not found
        ValueError: If barcode matching fails, ID parsing fails, or neither label_barcode nor seq_id provided
    """
    if not label_barcode and not seq_id:
        raise ValueError("Either label_barcode or seq_id must be provided")

    # Load metadata files
    seq_ids_path = f"{base_maple_dir}/mutation_data/{tag}/{barcode}/{tag}_{barcode}_seq-IDs.csv"
    genotypes_path = f"{base_maple_dir}/mutation_data/{tag}/{barcode}/{tag}_{barcode}_genotypes.csv"

    if not os.path.exists(seq_ids_path):
        raise FileNotFoundError(f"Metadata not found at {seq_ids_path}")

    seq_IDs = pd.read_csv(seq_ids_path)
    genotypes = pd.read_csv(genotypes_path)
    mut_data_df = pd.merge(seq_IDs, genotypes, on=['genotype_ID', 'reference_name'], how='left')

    # Identify target consensus sequence
    if seq_id:
        # Use seq_id directly
        matches = mut_data_df[mut_data_df['seq_ID'] == seq_id]
        if len(matches) != 1:
            raise ValueError(f"Expected 1 match for seq_ID '{seq_id}', found {len(matches)}")
        full_id = seq_id
    else:
        # Look up seq_id using label_barcode
        matches = mut_data_df[mut_data_df['barcode(s)'] == label_barcode]
        if len(matches) != 1:
            raise ValueError(f"Expected 1 match for barcode '{label_barcode}', found {len(matches)}")
        full_id = matches.iloc[0]['seq_ID']

    consensus_id_clean = full_id[:-2] if full_id.endswith('_0') else full_id

    try:
        parts = consensus_id_clean.split('|ref-')
        if len(parts) == 2:
            umi_tag, reference_name = parts[0], parts[1]
        else:
            parts = consensus_id_clean.split('|')
            umi_tag = parts[0]
            reference_name = parts[1].replace('ref-', '')
        target_unique_id = int(umi_tag.split('-')[1])
    except Exception as e:
        raise ValueError(f"Could not parse ID {consensus_id_clean}: {e}")

    # Load UMI groups log to find original reads
    groups_log_filename = f"{tag}_UMI-groups-log.csv"
    groups_log_paths = [
        f"{base_maple_dir}/sequences/UMI/{groups_log_filename}",
        f"{base_maple_dir}/{groups_log_filename}",
        groups_log_filename
    ]

    target_read_ids = set()
    for path in groups_log_paths:
        if os.path.exists(path):
            print(f"  Loading groups log from {path}")
            groups_df = pd.read_csv(path)
            subset = groups_df[
                (groups_df['unique_id'] == target_unique_id) &
                (groups_df['reference_name'] == reference_name)
            ]
            target_read_ids = set(subset['read_id'].astype(str))
            print(f"  Found {len(target_read_ids)} original reads for {umi_tag}")
            break

    if not target_read_ids:
        print(f"  Warning: No read IDs found in groups log, will use prefix matching")

    # Load reference sequences
    tags_df = pd.read_csv(f"{base_maple_dir}/metadata/tags.csv")
    ref_csv_name = tags_df.loc[tags_df['tag'] == tag, 'reference_csv'].values[0]
    ref_df = pd.read_csv(f"{base_maple_dir}/metadata/{ref_csv_name}", index_col='ref_seq_name')

    ref_seq_alignment = ref_df.loc[reference_name]['ref_seq_alignment'].upper()

    # Try getting separate coding sequence, fallback to alignment if missing
    try:
        ref_seq_coding = ref_df.loc[reference_name]['ref_seq_coding'].upper()
    except KeyError:
        print("  Warning: 'ref_seq_coding' not found, using alignment sequence")
        ref_seq_coding = ref_seq_alignment.replace("-", "")

    # Extract reads from BAM file
    bam_file = f"{base_maple_dir}/sequences/UMI/{tag}_pre-consensus.bam"
    aligned_data = []
    target_prefix = f"{umi_tag}|ref-{reference_name}"

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (target_read_ids and read.query_name in target_read_ids) or \
               (not target_read_ids and read.query_name.startswith(target_prefix)):
                seq_str, ins_locs = reconstruct_alignment_from_bam(read, len(ref_seq_alignment))
                aligned_data.append((seq_str, ins_locs, read.query_name))

    # Prepare metadata
    meta = {
        "UMI_ID": umi_tag,
        "Ref": reference_name,
        "Nonsyn": matches.iloc[0].get('AA_substitutions_nonsynonymous', 'None'),
        "Syn": matches.iloc[0].get('AA_substitutions_synonymous', 'None')
    }

    return aligned_data, ref_seq_alignment, ref_seq_coding, meta


def generate_alignment_plot(data_tuple, output_img, output_fasta):
    """
    Generates alignment visualization and saves FASTA file.

    Creates a comprehensive alignment plot showing:
    - Protein translation track with amino acid indices
    - Reference sequence
    - Up to 30 aligned reads with variants highlighted

    Args:
        data_tuple: Tuple of (reads_data, ref_seq_aln, ref_seq_coding, metadata)
        output_img: Path for output PNG image
        output_fasta: Path for output FASTA file
    """
    reads_data, ref_seq_aln, ref_seq_coding, meta = data_tuple

    if not reads_data:
        print(f"  Warning: No reads found for {meta['UMI_ID']}")
        return

    reads_to_plot = reads_data[:30]

    # Save FASTA file
    with open(output_fasta, 'w') as f:
        f.write(f">Reference\n{ref_seq_aln}\n")
        for seq_str, _, name in reads_data:
            f.write(f">{name}\n{seq_str}\n")
    print(f"  Saved FASTA: {output_fasta}")

    # Setup plot with dynamic height and fixed width
    fig_height = 0.5 * len(reads_to_plot) + 3
    fig, ax = plt.subplots(figsize=(40, fig_height))

    ax.set_title(
        f"Consensus: {meta['UMI_ID']} | Ref: {meta['Ref']}\n"
        f"Exp. Nonsyn: {meta['Nonsyn']} | Syn: {meta['Syn']}",
        fontsize=16, pad=30
    )

    FONT_SIZE = 10

    # Plot protein track at top
    prot_mapping = calculate_protein_mapping(ref_seq_aln, ref_seq_coding)
    prot_y_pos = len(reads_to_plot) + 1.2

    for item in prot_mapping:
        x_pos = item['plot_idx']
        aa_char = item['aa']
        aa_idx = item['aa_idx']

        # Center amino acid letter over codon
        ax.text(
            x_pos + 1, prot_y_pos, str(aa_char),
            ha='center', va='bottom', color='darkred',
            fontweight='bold', fontsize=FONT_SIZE + 2
        )

        # Show index vertically to save space
        ax.text(
            x_pos + 1, prot_y_pos + 0.5, str(aa_idx),
            ha='center', va='bottom', color='grey',
            fontsize=8, rotation=90
        )

    # Plot reference nucleotides
    ref_y_pos = len(reads_to_plot)
    for i, char in enumerate(ref_seq_aln):
        ax.text(
            i, ref_y_pos, char,
            ha='center', va='center',
            fontweight='bold', fontsize=FONT_SIZE, color='blue'
        )

    # Plot aligned reads
    for row, (seq_str, ins_locs, _) in enumerate(reads_to_plot):
        y = len(reads_to_plot) - 1 - row
        for col, char in enumerate(seq_str):
            # Mark insertions
            if col in ins_locs:
                ax.text(
                    col + 0.5, y, "^",
                    ha='center', va='center',
                    color='magenta', fontsize=14, fontweight='bold', zorder=10
                )

            if char == ' ':
                continue

            # Highlight deletions and mismatches
            if char == '-':
                rect = patches.Rectangle(
                    (col - 0.5, y - 0.5), 1, 1,
                    facecolor='#ffcccc', edgecolor='none'
                )
                ax.add_patch(rect)
            elif char != ref_seq_aln[col]:
                rect = patches.Rectangle(
                    (col - 0.5, y - 0.5), 1, 1,
                    facecolor='#ffffcc', edgecolor='none'
                )
                ax.add_patch(rect)
                ax.text(
                    col, y, char,
                    ha='center', va='center',
                    color='black', fontweight='bold', fontsize=FONT_SIZE
                )
            else:
                ax.text(
                    col, y, ".",
                    ha='center', va='center',
                    color='lightgrey', fontsize=FONT_SIZE
                )

    # Format plot
    ax.set_xlim(-0.5, len(ref_seq_aln) - 0.5)
    ax.set_ylim(-1, len(reads_to_plot) + 3)
    ax.set_yticks(range(len(reads_to_plot)))

    # Show full read IDs on y-axis (last part after splitting by '|')
    read_ids = [n.split('|')[-1] for _, _, n in reads_to_plot]
    ax.set_yticklabels(read_ids[::-1], fontsize=8)
    ax.set_xlabel("Reference Coordinate (bp)")

    plt.tight_layout()
    plt.savefig(output_img, dpi=150)
    plt.close()
    print(f"  Saved plot: {output_img}")


def process_csv_file(csv_path, input_dir, output_dir):
    """
    Processes multiple samples from a CSV file.

    Args:
        csv_path: Path to CSV file with columns: tag, barcode_group, and either
                  label_barcode or seq_ID (seq_ID takes precedence if both present)
        input_dir: Directory containing MAPLE pipeline outputs
        output_dir: Directory for output files
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read CSV file
    try:
        samples_df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        sys.exit(1)

    # Validate required columns
    required_cols = ['tag', 'barcode_group']
    missing_cols = [col for col in required_cols if col not in samples_df.columns]
    if missing_cols:
        print(f"Error: CSV missing required columns: {missing_cols}")
        print(f"Required columns: {required_cols}")
        sys.exit(1)

    # Check if either seq_ID or label_barcode is present
    has_seq_id = 'seq_ID' in samples_df.columns
    has_label_barcode = 'label_barcode' in samples_df.columns

    if not has_seq_id and not has_label_barcode:
        print("Error: CSV must have either 'seq_ID' or 'label_barcode' column")
        sys.exit(1)

    print(f"Processing {len(samples_df)} samples from {csv_path}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Using: {'seq_ID' if has_seq_id else 'label_barcode'}")
    print()

    # Process each sample
    success_count = 0
    error_count = 0

    for idx, row in samples_df.iterrows():
        tag = row['tag']
        barcode_group = row['barcode_group']

        # Get seq_id or label_barcode
        seq_id = row.get('seq_ID', None) if has_seq_id else None
        label_barcode = row.get('label_barcode', None) if has_label_barcode else None

        # Determine what to display and use for output filename
        if seq_id:
            identifier = seq_id
            output_identifier = seq_id.replace('|', '_').replace('/', '_')
            print(f"[{idx + 1}/{len(samples_df)}] Processing: {tag} / {barcode_group} / seq_ID: {seq_id}")
        else:
            identifier = label_barcode
            output_identifier = label_barcode
            print(f"[{idx + 1}/{len(samples_df)}] Processing: {tag} / {barcode_group} / {label_barcode}")

        try:
            # Get aligned sequences
            data = get_original_sequences_aligned(
                tag, barcode_group, input_dir,
                label_barcode=label_barcode,
                seq_id=seq_id
            )

            # Generate output filenames based on CSV inputs
            output_base = f"{tag}_{barcode_group}_{output_identifier}"
            output_fasta = os.path.join(output_dir, f"{output_base}.fasta")
            output_png = os.path.join(output_dir, f"{output_base}.png")

            # Generate plot and save files
            generate_alignment_plot(data, output_png, output_fasta)

            success_count += 1
            print()

        except Exception as e:
            print(f"  ERROR: {e}")
            error_count += 1
            print()

    print(f"Completed: {success_count} successful, {error_count} errors")


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Generate consensus read alignment visualizations from CSV input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example CSV format (using label_barcode):
  tag,barcode_group,label_barcode
  BLF1531-R1,A,B5
  BLF1531-R1,B,C3

Example CSV format (using seq_ID directly):
  tag,barcode_group,seq_ID
  BLF1531-R1,A,UMI-12345|ref-gene1_0
  BLF1531-R1,B,UMI-67890|ref-gene2_0

Example usage:
  python compile_original_reads.py samples.csv
  python compile_original_reads.py samples.csv -d /path/to/maple/output
  python compile_original_reads.py samples.csv -o results/
        """
    )

    parser.add_argument(
        'csv_file',
        help='CSV file with columns: tag, barcode_group, and either label_barcode or seq_ID'
    )

    parser.add_argument(
        '-d', '--directory',
        default='.',
        help='Directory containing MAPLE pipeline outputs (default: current directory)'
    )

    parser.add_argument(
        '-o', '--output',
        default='outputs',
        help='Output directory for FASTA and PNG files (default: outputs)'
    )

    args = parser.parse_args()

    # Validate input files exist
    if not os.path.exists(args.csv_file):
        print(f"Error: CSV file not found: {args.csv_file}")
        sys.exit(1)

    if not os.path.exists(args.directory):
        print(f"Error: Input directory not found: {args.directory}")
        sys.exit(1)

    # Process the CSV file
    process_csv_file(args.csv_file, args.directory, args.output)


if __name__ == '__main__':
    main()
