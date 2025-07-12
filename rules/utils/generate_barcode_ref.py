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
        dict: {barcode_type: {sequence: count}}
    """
    # Load reference sequence
    reference = list(SeqIO.parse(reference_fasta, 'fasta'))[0]
    reference_seq = str(reference.seq).upper()
    
    # Initialize barcode counts for types that need generation
    barcode_counts = {}
    for barcode_type, info in barcode_info.items():
        if info.get('generate', False) and not os.path.exists(info['fasta']):
            barcode_counts[barcode_type] = Counter()
    
    if not barcode_counts:
        return {}
    
    # Open BAM file and process reads
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    
    for bam_entry in bam_file.fetch(reference.id):
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
    
    bam_file.close()
    
    return barcode_counts


def filter_barcodes_by_hamming_distance(barcode_df, min_hamming_distance):
    """
    Filter barcodes to ensure minimum hamming distance between selected sequences.
    
    Args:
        barcode_df: DataFrame with columns ['barcode', 'count'] sorted by count
        min_hamming_distance: Minimum hamming distance required
        
    Returns:
        list: Selected barcode sequences
    """
    selected_barcodes = []
    
    for _, row in barcode_df.iterrows():
        barcode = row['barcode']
        
        # Check if this barcode is far enough from all selected barcodes
        too_close = False
        for selected in selected_barcodes:
            if hamming_distance(barcode, selected) <= min_hamming_distance:
                too_close = True
                break
        
        if not too_close:
            selected_barcodes.append(barcode)
    
    return selected_barcodes


def write_barcode_fasta(barcode_counts, barcode_info, output_flag):
    """
    Write barcode FASTA files for each barcode type.
    
    Args:
        barcode_counts: Dictionary of barcode counts {type: {seq: count}}
        barcode_info: Dictionary of barcode information
        output_flag: Path to output flag file
        
    Returns:
        None
    """
    if not barcode_counts:
        # No barcodes to generate, just create flag file
        Path(output_flag).touch()
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
            selected_barcodes = filter_barcodes_by_hamming_distance(df, hamming_dist)
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
    
    # Create flag file
    Path(output_flag).touch()


def main():
    # Parse arguments
    if 'snakemake' in globals():
        # Running as snakemake script
        args = type('Args', (), {
            'bam': snakemake.input.bam,
            'reference': snakemake.params.reference,
            'barcode_info': snakemake.input.barcode_info,
            'tag': snakemake.wildcards.tag,
            'output_flag': snakemake.output.flag
        })()
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--bam', required=True, help='Input BAM file')
        parser.add_argument('--reference', required=True, help='Reference FASTA file')
        parser.add_argument('--barcode_info', required=True, help='Barcode info CSV file')
        parser.add_argument('--tag', required=True, help='Sample tag')
        parser.add_argument('--output_flag', required=True, help='Output flag file')
        args = parser.parse_args()
    
    # Load barcode information
    barcode_info = load_csv_as_dict(
        args.barcode_info,
        required=['barcode_name', 'context', 'fasta'],
        tag=args.tag,
        defaults={'reverse_complement': False, 'generate': False, 'hamming_distance': 0}
    )
    
    # Validate barcode info
    metadata_dir = Path(args.barcode_info).parts[0]
    for bc in barcode_info:
        barcode_info, _ = validate_bc(barcode_info, bc, metadata_dir=metadata_dir)
    
    # Find barcodes in BAM file
    barcode_counts = find_barcodes_in_bam(args.bam, args.reference, barcode_info)
    
    # Write barcode FASTA files
    write_barcode_fasta(barcode_counts, barcode_info, args.output_flag)


if __name__ == '__main__':
    main()