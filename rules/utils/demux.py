#!/usr/bin/env python3
"""
Demultiplexing of BAM files based on barcode sequences.

This script processes aligned BAM files to demultiplex reads based on one or more barcode types.
Each barcode type is defined by:
- A sequence context (e.g., "ATCGNNNNCCGA") where N marks the barcode position
- A FASTA file containing valid barcode sequences
- Optional parameters like reverse complement and hamming distance tolerance

Key features:
1. Multiple barcode types: Supports any number of different barcodes per read
2. Partition barcode groups: Barcode combinations that determine output file partitioning
3. Label barcode groups: Barcode combinations that add labels to sequences without partitioning
4. Barcode modes: Each barcode type can be used for partitioning, labeling, or both
5. Hamming distance tolerance: Allows imperfect barcode matches within specified distance
6. IUPAC ambiguity codes: Treats N,Y,R,S,W,K,M,B,D,H,V as wildcards in contexts
7. Quality filtering: Removes files with too few reads or failed barcodes
8. Group pooling: Multiple barcode combinations can map to the same output file

The script can run standalone or as part of a Snakemake pipeline.
"""

import os
from pathlib import Path
import sys
import argparse
import itertools
from collections import Counter
from common import load_csv_as_dict, validate_bc

import numpy as np
import pandas as pd
import pysam
from Bio import Seq, SeqIO
from sklearn.metrics import pairwise_distances


# ============================================================================
# Common utility functions (can be moved to common.py)
# ============================================================================

def hamming_distance(string1, string2):
    """
    Calculate hamming distance between two strings.
    
    Args:
        string1: First string
        string2: Second string
        
    Returns:
        int: Number of differing positions
    """
    return sum(c1 != c2 for c1, c2 in zip(string1.upper(), string2.upper()))


def create_barcodes_dict(barcode_fasta, reverse_complement=False):
    """
    Create a dictionary of barcodes from fasta file.
    
    Args:
        barcode_fasta: Path to barcode fasta file
        reverse_complement: If True, store reverse complements of sequences
        
    Returns:
        dict: {sequence: barcode_name}
    """
    barcode_dict = {}
    for entry in SeqIO.parse(barcode_fasta, 'fasta'):
        barcode_name = entry.id
        sequence = str(entry.seq).upper()
        if reverse_complement:
            sequence = str(Seq.reverse_complement(sequence))
        barcode_dict[sequence] = barcode_name
    return barcode_dict


def encode_sequence(seq):
    """Encode DNA sequence as integers for hamming distance computation."""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3, '-': 4}
    return np.array([mapping[c] for c in seq])


def find_barcode_hamming(query_seq, barcode_matrix, barcode_names, max_distance):
    """
    Find closest barcode within hamming distance using sklearn.

    Args:
        query_seq: Query barcode sequence
        barcode_matrix: Pre-vectorized N×L integer array of reference barcodes
        barcode_names: List of barcode names corresponding to matrix rows
        max_distance: Maximum hamming distance allowed

    Returns:
        str or None: Barcode name if match found within distance, None otherwise
    """
    query_array = encode_sequence(query_seq).reshape(1, -1)
    distances = pairwise_distances(query_array, barcode_matrix, metric='hamming')
    distances_int = (distances[0] * len(query_seq)).astype(int)

    min_idx = np.argmin(distances_int)
    return barcode_names[min_idx] if distances_int[min_idx] <= max_distance else None


# ============================================================================
# Core demultiplexing functions
# ============================================================================

def load_csv_with_tag_filter(csv_path, key_column, tag=None):
    """
    Load CSV file with optional tag filtering, preserving duplicate keys.
    
    Args:
        csv_path: Path to CSV file
        key_column: Column name to use as group identifier
        tag: Optional tag to filter by
        
    Returns:
        list: List of (key, row_dict) tuples to preserve duplicates
    """
    if not csv_path or not os.path.exists(csv_path):
        return []
    
    df = pd.read_csv(csv_path)
    
    # Filter by tag if specified and tag column exists
    if tag and 'tag' in df.columns:
        df = df[df['tag'] == tag]
    
    result = []
    for _, row in df.iterrows():
        key = row[key_column]
        row_dict = row.to_dict()
        row_dict.pop(key_column, None)
        row_dict.pop('tag', None)
        # Remove any NaN values
        row_dict = {k: v for k, v in row_dict.items() if pd.notna(v)}
        if row_dict:
            result.append((key, row_dict))
    
    return result


def load_barcode_info(barcode_info_csv, tag=None):
    """
    Load barcode information from CSV file and validate it is properly defined.
    
    Args:
        barcode_info_csv: Path to CSV file
        tag: Optional tag to filter by
        
    Returns:
        dict: {barcode_type: info_dict}
    """

    # parent folder of barcode_info_csv is the metadata dir
    metadata_dir = Path(barcode_info_csv).parts[0]

    barcode_info = load_csv_as_dict(
        barcode_info_csv,
        required=['barcode_name', 'context', 'fasta', 'reverse_complement'],
        tag=tag,
        defaults={'reverse_complement': False, 'label_only': False, 'generate': False, 'hamming_distance': 0}
    )
    
    # Validate each barcode type
    for bc in barcode_info:
        barcode_info, _ = validate_bc(barcode_info, bc, metadata_dir=metadata_dir)
    
    return barcode_info


def validate_barcodes(barcode_info, reference_seq):
    """
    Validate that all barcode contexts can be found in reference.
    
    Args:
        barcode_info: Dictionary of barcode information
        reference_seq: Reference sequence string
        
    Raises:
        ValueError: If context not found in reference
    """
    iupac_codes = 'NYRSWKMBDHV'
    
    for barcode_type, info in barcode_info.items():
        context = info['context'].upper()
        
        # Convert IUPAC codes to regex pattern
        import re
        pattern = ''
        for char in context:
            if char in iupac_codes:
                pattern += '.'
            else:
                pattern += char
        
        if not re.search(pattern, reference_seq):
            raise ValueError(
                f"Barcode context not found in reference sequence.\n"
                f"Barcode type: {barcode_type}\n"
                f"Context: {context}"
            )


def build_barcode_dicts(barcode_info):
    """
    Build dictionaries for barcode lookup.

    Args:
        barcode_info: Dictionary of barcode information

    Returns:
        tuple: (barcode_dicts, hamming_lookups)
    """
    barcode_dicts = {}
    hamming_lookups = {}

    for barcode_type, info in barcode_info.items():
        fasta_path = info['fasta']
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"Barcode fasta not found: {fasta_path}")

        barcode_dict = create_barcodes_dict(
            fasta_path,
            info.get('reverse_complement', False)
        )
        barcode_dicts[barcode_type] = barcode_dict

        # Cache vectorized sequences for hamming distance
        hamming_dist = info.get('hamming_distance', 0)
        if hamming_dist > 0:
            sequences = list(barcode_dict.keys())
            names = list(barcode_dict.values())
            # Pre-vectorize into N×L integer array - CACHED!
            barcode_matrix = np.array([encode_sequence(s) for s in sequences])
            hamming_lookups[barcode_type] = {
                'matrix': barcode_matrix,
                'names': names,
                'max_distance': hamming_dist,
                'cached_matches': {}  # Separate dict for hamming-found matches
            }

    return barcode_dicts, hamming_lookups


def build_reference_alignment(bam_entry, reference_seq):
    """
    Build reference alignment string with indels.
    
    Args:
        bam_entry: pysam AlignmentFile entry
        reference_seq: Reference sequence string
        
    Returns:
        str: Aligned reference sequence with gaps
    """
    index = bam_entry.reference_start
    ref_alignment = ''
    
    for op, length in bam_entry.cigartuples:
        if op == 0:  # Match
            ref_alignment += reference_seq[index:index + length]
            index += length
        elif op == 1:  # Insertion
            ref_alignment += '-' * length
        elif op == 2:  # Deletion
            index += length
    
    return ref_alignment


def find_barcode_position_in_alignment(ref_alignment, context):
    """
    Find barcode position in aligned sequence.
    
    Args:
        ref_alignment: Aligned reference sequence
        context: Barcode context with Ns marking position
        
    Returns:
        tuple: (start_pos, end_pos) or (None, error_reason)
    """
    # Simple string find (matching original logic)
    location = ref_alignment.find(context)
    if location == -1:
        return None, 'context_not_present_in_reference_sequence'
    
    n_start = context.find('N')
    n_end = context.rfind('N') + 1
    
    barcode_start = location + n_start
    barcode_end = location + n_end
    
    # Check for indels near barcode (matching original logic)
    if ((barcode_start > 0 and ref_alignment[barcode_start-1] == '-') or 
        (barcode_end < len(ref_alignment) and ref_alignment[barcode_end] == '-')):
        return None, 'low_confidence_barcode_identification'
    
    return barcode_start, barcode_end


def identify_barcodes(bam_entry, reference_seq, barcode_info, barcode_dicts, hamming_lookups):
    """
    Identify all barcodes in a sequence.

    Args:
        bam_entry: pysam AlignmentFile entry
        reference_seq: Reference sequence string
        barcode_info: Dictionary of barcode information
        barcode_dicts: Dictionary of barcode lookup dicts
        hamming_lookups: Dictionary of hamming lookup info with cached matrices

    Returns:
        tuple: (sequence_barcodes_dict, barcode_names_list, stats_array)
    """
    ref_alignment = build_reference_alignment(bam_entry, reference_seq)

    sequence_barcodes = {}
    barcode_names = []
    stats_list = []

    for barcode_type, info in barcode_info.items():
        context = info['context'].upper()
        barcode_dict = barcode_dicts[barcode_type]
        hamming_lookup = hamming_lookups.get(barcode_type)

        start, end = find_barcode_position_in_alignment(ref_alignment, context)

        # Initialize stats
        not_exact_match = 0
        failure_reasons = {
            'context_not_present_in_reference_sequence': 0,
            'barcode_not_in_fasta': 0,
            'low_confidence_barcode_identification': 0
        }

        if isinstance(start, int):
            # Extract barcode from query sequence
            barcode_seq = bam_entry.query_alignment_sequence[start:end]

            # Try exact match, then cached hamming matches, then compute hamming
            barcode_name = barcode_dict.get(barcode_seq)
            if barcode_name is None and hamming_lookup:
                barcode_name = hamming_lookup['cached_matches'].get(barcode_seq)
                if barcode_name:
                    not_exact_match = 1
                else:
                    barcode_name = find_barcode_hamming(
                        barcode_seq, hamming_lookup['matrix'],
                        hamming_lookup['names'], hamming_lookup['max_distance']
                    )
                    if barcode_name:
                        not_exact_match = 1
                        hamming_lookup['cached_matches'][barcode_seq] = barcode_name
            if barcode_name is None:
                barcode_name = 'fail'
                failure_reasons['barcode_not_in_fasta'] = 1
        else:
            # start is None, end contains error reason
            barcode_name = 'fail'
            error_key = 'context_not_present_in_reference_sequence' if end == 'context_not_present' else 'low_confidence_barcode_identification'
            failure_reasons[error_key] = 1

        sequence_barcodes[barcode_type] = barcode_name
        barcode_names.append(barcode_name)
        # Only append stats, not barcode names
        stats_list.extend([not_exact_match] + [failure_reasons['context_not_present_in_reference_sequence'],
                                                failure_reasons['barcode_not_in_fasta'],
                                                failure_reasons['low_confidence_barcode_identification']])

    return sequence_barcodes, barcode_names, np.array(stats_list)


def get_group_for_barcodes(sequence_barcodes, groups_list, barcode_types):
    """
    Find matching group for barcode combination.
    
    Args:
        sequence_barcodes: Dict of identified barcodes {type: name}
        groups_list: List of (group_name, group_def) tuples
        barcode_types: List of barcode types to check
        
    Returns:
        tuple: (group_name or None, remaining_barcode_names)
    """
    # Get barcodes for the specified types
    group_barcodes = {bt: sequence_barcodes.get(bt, 'fail') for bt in barcode_types}
    remaining_barcodes = {bt: name for bt, name in sequence_barcodes.items() 
                         if bt not in barcode_types}
    
    # Check if any barcode failed
    if 'fail' in group_barcodes.values() and barcode_types:
        return None, list(sequence_barcodes.values())
    
    # Look for matching group
    for group_name, group_def in groups_list:
        if all(group_def.get(bt) == group_barcodes.get(bt) for bt in barcode_types if bt in group_def):
            return group_name, list(remaining_barcodes.values())
    
    return None, list(sequence_barcodes.values())


def demultiplex_bam(bam_path, output_dir, barcode_info, partition_groups, 
                    label_groups, tag, reference_fasta, 
                    demux_threshold, demux_screen_failures, demux_screen_no_group):
    """
    Main demultiplexing function.
    
    Args:
        bam_path: Input BAM file path
        output_dir: Output directory path
        barcode_info: Barcode information dictionary
        partition_groups: List of partition group tuples
        label_groups: List of label group tuples
        tag: Sample tag
        reference_fasta: Reference FASTA file path
        demux_threshold: Minimum threshold for demultiplexing. If greater than 1, it is treated as an absolute count.
        demux_screen_failures: Whether to screen failed barcodes
        demux_screen_no_group: Whether to screen ungrouped files
        
    Returns:
        None (writes output files and statistics)
    """
    # Load reference sequence
    reference = list(SeqIO.parse(reference_fasta, 'fasta'))[0]
    reference_seq = str(reference.seq).upper()
    
    # Validate barcodes
    validate_barcodes(barcode_info, reference_seq)
    
    # Build barcode dictionaries
    barcode_dicts, hamming_lookups = build_barcode_dicts(barcode_info)
    
    # Categorize barcode types
    partition_bc_types = []
    label_bc_types = []
    
    # Identify which barcodes are used in groups
    partition_group_barcodes = set()
    for _, group_def in partition_groups:
        partition_group_barcodes.update(group_def.keys())
    
    label_group_barcodes = set()
    for _, group_def in label_groups:
        label_group_barcodes.update(group_def.keys())
    
    # Categorize each barcode type
    for barcode_type in barcode_info:
        if barcode_info[barcode_type].get('label_only', False):
            label_bc_types.append(barcode_type)
        else:
            partition_bc_types.append(barcode_type)
    
    # Open BAM file
    bam_file = pysam.AlignmentFile(bam_path, 'rb')
    output_files = {}
    row_counts = Counter()
    
    # Process each read
    for bam_entry in bam_file.fetch(reference.id):
        # Identify barcodes
        sequence_barcodes, barcode_names, stats_array = identify_barcodes(
            bam_entry, reference_seq, barcode_info, barcode_dicts, hamming_lookups
        )
        
        # Handle label barcodes - add as BAM tag
        label_barcodes = {bt: sequence_barcodes[bt] for bt in label_bc_types 
                         if bt in sequence_barcodes}
        
        label_output = ''  # Track label combination for stats
        if label_barcodes:
            # Find matching label group
            label_name, remaining = get_group_for_barcodes(
                label_barcodes, label_groups, 
                [bt for bt in label_bc_types if bt in label_group_barcodes]
            )
            
            if label_name and remaining:
                label_tag = f"{label_name}_{'_'.join(remaining)}"
                label_output = f"{label_name}-{'-'.join(remaining)}"
            elif label_name:
                label_tag = label_name
                label_output = label_name
            else:
                label_tag = '_'.join(label_barcodes.values())
                label_output = '-'.join(label_barcodes.values())
            
            bam_entry.set_tag('BC', label_tag)
        else:
            label_output = 'none'  # No label barcodes or all failed
        
        # Get output filename using partition barcodes
        partition_barcodes = {bt: sequence_barcodes[bt] for bt in partition_bc_types 
                             if bt in sequence_barcodes}
        
        group_name, remaining = get_group_for_barcodes(
            partition_barcodes, partition_groups,
            [bt for bt in partition_bc_types if bt in partition_group_barcodes]
        )
        
        if group_name and remaining:
            output_prefix = f"{group_name}-{'-'.join(remaining)}"
        elif group_name:
            output_prefix = group_name
        else:
            output_prefix = '-'.join(partition_barcodes.values()) if partition_barcodes else 'all'
        
        is_grouped = group_name is not None
        
        # Create output file if needed
        if output_prefix not in output_files:
            output_path = os.path.join(output_dir, f'{tag}_{output_prefix}.bam')
            output_files[output_prefix] = pysam.AlignmentFile(
                output_path, 'wb', template=bam_file
            )
        
        # Write to output file
        output_files[output_prefix].write(bam_entry)
        
        # Update statistics - include label_output in the key
        # Order: tag, output_prefix, is_grouped, ...barcode_names, label_output
        stats_key = tuple([tag, output_prefix, is_grouped] + barcode_names + [label_output])
        stats_array = np.insert(stats_array, 0, 1)  # Add count
        row_counts[stats_key] += stats_array
    
    # Close files
    bam_file.close()
    for f in output_files.values():
        f.close()
    
    # Generate and save statistics
    stats_df = generate_statistics_dataframe(row_counts, barcode_info, tag, partition_groups, label_groups)
    
    # Filter output files
    if barcode_info:
        filter_output_files(stats_df, output_dir, tag, barcode_info,
                        demux_threshold, demux_screen_failures, demux_screen_no_group)
    
    # Save statistics
    stats_df.drop('named_by_group', axis=1, inplace=True)
    stats_df.to_csv(os.path.join(output_dir, f'{tag}_demux-stats.csv'), index=False)


def generate_statistics_dataframe(row_counts, barcode_info, tag, partition_groups, label_groups):
    """
    Generate statistics dataframe from counting results.
    
    Args:
        row_counts: Counter object with statistics
        barcode_info: Barcode information dictionary
        tag: Sample tag
        partition_groups: List of partition group tuples
        label_groups: List of label group tuples
        
    Returns:
        pandas.DataFrame: Statistics dataframe
    """
    # Build column names - ensuring label_barcodes is right before barcodes_count
    columns = ['tag', 'output_file_barcodes', 'named_by_group']
    
    # Add barcode type columns
    for barcode_type in barcode_info:
        columns.append(barcode_type)
    
    # Add label_barcodes right before barcodes_count
    columns.append('label_barcodes')
    columns.append('barcodes_count')
    
    # Add failure columns for each barcode type
    for barcode_type in barcode_info:
        columns.extend([
            f'{barcode_type}:not_exact_match',
            f'{barcode_type}_failed:context_not_present_in_reference_sequence',
            f'{barcode_type}_failed:barcode_not_in_fasta',
            f'{barcode_type}_failed:low_confidence_barcode_identification'
        ])
    
    # Convert counter to rows
    rows = []
    for key, counts in row_counts.items():
        # key format: (tag, output_file_barcodes, named_by_group, ...barcode_names, label_barcodes)
        # counts format: [barcodes_count, ...failure stats for each barcode type...]
        
        tag_val, output_file, is_grouped = key[:3]
        num_barcodes = len(barcode_info)
        barcode_names = key[3:3+num_barcodes]  # Extract exact number of barcode names
        label_output = key[3+num_barcodes]     # Label output comes after barcode names
        
        row = [tag_val, output_file, is_grouped]
        row.extend(barcode_names)      # Individual barcode names
        row.append(label_output)        # label_barcodes
        row.extend(counts.tolist())     # barcodes_count + all failure stats
        rows.append(row)
    
    if not rows:
        # Create row for each expected partition group with 0 counts
        for group_name, _ in partition_groups:
            row = [tag, group_name, True]
            row.extend(['n/a'] * len(barcode_info))  # barcode names
            row.append('none')  # label_barcodes
            row.append(0)  # barcodes_count
            row.extend([0] * (4 * len(barcode_info)))  # failure stats
            rows.append(row)
        # If no groups, create default empty row
        if not rows:
            row = [tag, 'all', False]
            row.extend(['fail'] * len(barcode_info))  # barcode names
            row.append('none')  # label_barcodes
            row.append(0)  # barcodes_count
            row.extend([0] * (4 * len(barcode_info)))  # failure stats
            rows.append(row)
    
    df = pd.DataFrame(rows, columns=columns)
    
    # Add file counts (demuxed_count)
    file_counts = df.groupby('output_file_barcodes')['barcodes_count'].sum()
    df['demuxed_count'] = df['output_file_barcodes'].map(file_counts)
    
    # Build aggregation dict
    group_by_columns = ['tag', 'output_file_barcodes', 'named_by_group', 'demuxed_count', 'label_barcodes']
    group_by_columns.extend(barcode_info.keys())  # Add barcode type columns
    
    sum_columns = {'barcodes_count': 'sum'}
    for barcode_type in barcode_info:
        for suffix in [':not_exact_match', '_failed:context_not_present_in_reference_sequence',
                      '_failed:barcode_not_in_fasta', '_failed:low_confidence_barcode_identification']:
            sum_columns[f'{barcode_type}{suffix}'] = 'sum'
    
    # Group and aggregate
    df = df.groupby(group_by_columns).agg(sum_columns).reset_index()
    
    # Add rows for expected partition groups that have no sequences
    existing_groups = set(df[df['named_by_group'] == True]['output_file_barcodes'])
    expected_groups = {group_name for group_name, _ in partition_groups}
    missing_groups = expected_groups - existing_groups
    
    for group_name in missing_groups:
        new_row = {
            'tag': tag,
            'output_file_barcodes': group_name,
            'named_by_group': True,
            'label_barcodes': 'none',
            'demuxed_count': 0,
            'barcodes_count': 0
        }
        # Add barcode columns with 'n/a' values
        for barcode_type in barcode_info:
            new_row[barcode_type] = 'n/a'
            for suffix in [':not_exact_match', '_failed:context_not_present_in_reference_sequence',
                          '_failed:barcode_not_in_fasta', '_failed:low_confidence_barcode_identification']:
                new_row[f'{barcode_type}{suffix}'] = 0
        
        df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    
    # Add rows for expected label groups within each output file group
    for output_file in df['output_file_barcodes'].unique():
        output_df = df[df['output_file_barcodes'] == output_file]
        existing_labels = set(output_df['label_barcodes'])
        expected_labels = {group_name for group_name, _ in label_groups}
        
        # Get the demuxed_count for this output file
        output_demuxed_count = output_df['demuxed_count'].iloc[0] if len(output_df) > 0 else 0
        
        # Add 'none' if there are any rows without label groups
        if expected_labels and 'none' not in existing_labels and any(lb not in expected_labels for lb in existing_labels if lb != 'none'):
            expected_labels.add('none')
        
        missing_labels = expected_labels - existing_labels
        
        for label_name in missing_labels:
            # Get barcode values from an existing row for this output file
            sample_row = output_df.iloc[0] if len(output_df) > 0 else None
            
            new_row = {
                'tag': tag,
                'output_file_barcodes': output_file,
                'named_by_group': sample_row['named_by_group'] if sample_row is not None else False,
                'label_barcodes': label_name,
                'demuxed_count': output_demuxed_count,
                'barcodes_count': 0
            }
            
            # Add barcode columns - use same values as other rows for this output file, but 'n/a' for label barcodes
            for barcode_type in barcode_info:
                if sample_row is not None and barcode_type in sample_row:
                    if barcode_info[barcode_type].get('label_only', False):
                        new_row[barcode_type] = 'n/a'
                    else:
                        new_row[barcode_type] = sample_row[barcode_type]
                else:
                    new_row[barcode_type] = 'n/a'
                    
                for suffix in [':not_exact_match', '_failed:context_not_present_in_reference_sequence',
                              '_failed:barcode_not_in_fasta', '_failed:low_confidence_barcode_identification']:
                    new_row[f'{barcode_type}{suffix}'] = 0
            
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
    
    df.sort_values(['demuxed_count', 'barcodes_count'], ascending=False, inplace=True)

    # reorder to put the label_barcodes column right before barcodes_count
    bc_count_index = df.columns.get_loc('barcodes_count')
    cols = list(df.columns)
    cols.insert(bc_count_index - 1, cols.pop(cols.index('label_barcodes')))
    df = df[cols]
    
    return df


def filter_output_files(stats_df, output_dir, tag, barcode_info, 
                       demux_threshold, demux_screen_failures, demux_screen_no_group):
    """
    Move files that don't meet criteria to subdirectory.
    
    Args:
        stats_df: Statistics dataframe
        output_dir: Output directory path
        tag: Sample tag
        barcode_info: Barcode information dictionary
        demux_threshold: Minimum proportion threshold
        demux_screen_failures: Whether to screen failed barcodes
        demux_screen_no_group: Whether to screen ungrouped files
        
    Returns:
        None (moves files)
    """
    banish_dir = os.path.join(output_dir, 'no_subsequent_analysis')
    os.makedirs(banish_dir, exist_ok=True)
    
    total_sequences = stats_df['barcodes_count'].sum()
    
    # Determine barcode types
    partition_bc_types = [bt for bt, info in barcode_info.items() 
                         if not info.get('label_only', False)]
    label_bc_types = [bt for bt, info in barcode_info.items() 
                     if info.get('label_only', False)]
    
    banished = set()
    
    for _, row in stats_df.iterrows():
        if row['demuxed_count'] == 0 or row['output_file_barcodes'] in banished:
            continue
        
        should_banish = False
        
        # Check threshold
        if demux_threshold > 1: 
            threshold_absolute = demux_threshold
        else:
            threshold_absolute = demux_threshold * total_sequences
        if row['demuxed_count'] < threshold_absolute:
            should_banish = True
        
        # Check for failures - ANY failure in partition barcodes should banish
        if demux_screen_failures:
            # Check partition barcodes - any failure banishes the file
            for barcode_type in partition_bc_types:
                if barcode_type in row and row[barcode_type] == 'fail':
                    should_banish = True
                    break
            
            # Check label barcodes (only if all sequences have this barcode failed)
            for barcode_type in label_bc_types:
                if (barcode_type in row and row[barcode_type] == 'fail' and 
                    row['barcodes_count'] == row['demuxed_count']):
                    should_banish = True
                    break
        
        # Check group assignment
        if demux_screen_no_group and not row['named_by_group']:
            should_banish = True
        
        if should_banish:
            source = os.path.join(output_dir, f"{tag}_{row['output_file_barcodes']}.bam")
            dest = os.path.join(banish_dir, f"{tag}_{row['output_file_barcodes']}.bam")
            if os.path.exists(source):
                os.rename(source, dest)
            banished.add(row['output_file_barcodes'])


def main():
    # Parse arguments
    if 'snakemake' in globals():
        # Running as snakemake script
        args = type('Args', (), {
            'bam': snakemake.input.aln,
            'reference': snakemake.params.reference,
            'barcode_info': snakemake.input.barcode_info,
            'partition_barcode_groups': snakemake.input.get('partition_barcode_groups', None),
            'label_barcode_groups': snakemake.input.get('label_barcode_groups', None),
            'output_dir': os.path.dirname(str(snakemake.output.flag)),
            'stats': snakemake.output.stats,
            'tag': snakemake.wildcards.tag,
            'demux_threshold': snakemake.params.demux_threshold,
            'demux_screen_failures': snakemake.params.demux_screen_failures,
            'demux_screen_no_group': snakemake.params.demux_screen_no_group
        })()
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--bam', required=True, help='Input BAM file')
        parser.add_argument('--reference', required=True, help='Reference FASTA file')
        parser.add_argument('--barcode_info', required=True, help='Barcode info CSV file')
        parser.add_argument('--partition_barcode_groups', help='Partition barcode groups CSV file')
        parser.add_argument('--label_barcode_groups', help='Label barcode groups CSV file')
        parser.add_argument('--output_dir', required=True, help='Output directory')
        parser.add_argument('--stats', required=True, help='Output statistics file')
        parser.add_argument('--tag', required=True, help='Sample tag')
        parser.add_argument('--demux_threshold', type=float, default=0,
                          help='Minimum proportion of reads for file to be kept')
        parser.add_argument('--demux_screen_failures', action='store_true',
                          help='Remove files with failed barcodes')
        parser.add_argument('--demux_screen_no_group', action='store_true',
                          help='Remove files not assigned to a group')
        args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load inputs
    if args.barcode_info and os.path.exists(args.barcode_info):
        barcode_info = load_barcode_info(args.barcode_info, args.tag)
    else:
        barcode_info = {}
    partition_groups = load_csv_with_tag_filter(args.partition_barcode_groups, 'barcode_group', args.tag) if args.partition_barcode_groups else []
    label_groups = load_csv_with_tag_filter(args.label_barcode_groups, 'barcode_group', args.tag) if args.label_barcode_groups else []
    
    # Run demultiplexing
    demultiplex_bam(
        args.bam, args.output_dir, barcode_info, partition_groups,
        label_groups, args.tag, args.reference,
        args.demux_threshold, args.demux_screen_failures, args.demux_screen_no_group
    )


if __name__ == '__main__':
    main()