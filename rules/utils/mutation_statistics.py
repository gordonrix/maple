#!/usr/bin/env python3
"""
Mutation statistics calculation for maple pipeline.

This script processes mutation analysis output files to calculate comprehensive
statistics about mutations, including:
- Total and unique sequence counts
- NT and AA mutation rates, distributions, and spectra
- Transversion/transition ratios
- Mutation type frequencies (all 12 possible substitution types)
- Insertion and deletion statistics

The script can run standalone or as part of a Snakemake pipeline.
"""

import os
import sys
import argparse
import collections

import numpy as np
import pandas as pd
from Bio import SeqIO


def compute_mean_from_dist(dist):
    """Compute mean from pandas series distribution."""
    total = sum(n * count for n, count in enumerate(dist))
    return total / dist.sum() if total != 0 else 0


def compute_median_from_dist(dist):
    """Compute median from distribution."""
    if all(count == 0 for count in dist):
        return 0
    seq_list = [n for n, count in enumerate(dist) for _ in range(int(count))]
    return int(pd.Series(seq_list).median())


def mut_type(wt, muts):
    """
    Calculate mutation type counts for all 12 possible substitutions.

    Args:
        wt: Wild-type nucleotide (e.g., 'A1', 'C4')
        muts: Dict of mutation counts {'A': count, 'T': count, 'G': count, 'C': count}

    Returns:
        OrderedDict of mutation type counts
    """
    wt_nt = wt[0]
    nts = 'ATGC'
    out_dict = collections.OrderedDict()

    for nt in nts:
        for mut in nts:
            if mut != nt:
                out_dict[f'{nt}->{mut}'] = muts[mut] if (wt_nt in nts and nt == wt_nt) else 0
    return out_dict


def transversions_transitions(wt, muts):
    """Calculate transversions and transitions from mutation counts."""
    wt_nt = wt[0]

    if wt_nt in ['A', 'G']:
        return {'transversions': muts['T'] + muts['C'], 'transitions': muts['A'] + muts['G']}
    elif wt_nt in ['T', 'C']:
        return {'transversions': muts['A'] + muts['G'], 'transitions': muts['T'] + muts['C']}
    return {'transversions': 0, 'transitions': 0}


def build_input_file_dict(input_file_list):
    """
    Generate nested dictionary of input files organized by tag and barcode.

    Returns:
        dict: {tag: {barcode_group: {datatype: filename}}}
    """
    out_dict = {}
    for f in input_file_list:
        tag = f.split('_')[-3].split('/')[-1]
        barcodes = f.split('_')[-2]
        dtype = f.split('_')[-1].split('.')[0]

        out_dict.setdefault(tag, {}).setdefault(barcodes, {})[dtype] = f
    return out_dict


def calculate_statistics_for_sample(df_dict, tag, barcode_group,
                                     do_aa_analysis, mutations_frequencies_raw=False):
    """
    Calculate all statistics for a single tag/barcode combination.

    Args:
        df_dict: Dictionary of DataFrames with mutation data
        tag: Sample tag
        barcode_group: Barcode group name
        do_aa_analysis: Whether AA mutation analysis was performed
        mutations_frequencies_raw: Whether frequencies are raw proportions (default: False)
    """
    # Basic counts
    nt_dist = df_dict['NT-mutation-distribution']['total sequences']
    total_seqs = nt_dist.sum()
    fail_count = len(df_dict['failures'])
    reference_length = len(df_dict['NT-mutation-frequencies'])

    # Count unique genotypes
    num_genotypes = (df_dict['genotypes']
                     [df_dict['genotypes']['count'] > 0]
                     [['NT_substitutions', 'NT_insertions', 'NT_deletions']]
                     .drop_duplicates()
                     .shape[0])

    values_list = [tag, barcode_group, total_seqs, num_genotypes, fail_count]

    # AA mutation statistics
    if do_aa_analysis:
        aa_dist = df_dict['AA-mutation-distribution']['total sequences']
        if mutations_frequencies_raw:
            total_aa_muts = df_dict['AA-mutation-frequencies'].values.sum()
        else:
            total_aa_muts = int(round((df_dict['AA-mutation-frequencies'] * total_seqs).values.sum()))
        unique_aa_muts = int(df_dict['AA-mutation-frequencies'].where(
            df_dict['AA-mutation-frequencies'] == 0, 1).values.sum())
        values_list.extend([total_aa_muts, unique_aa_muts,
                           compute_mean_from_dist(aa_dist),
                           compute_median_from_dist(aa_dist)])
    else:
        values_list.extend(['N/A'] * 4)

    # NT mutation processing
    nt_muts = df_dict['NT-mutation-frequencies'].copy()
    if not mutations_frequencies_raw:
        nt_muts = np.rint(nt_muts * total_seqs).astype(int)

    nt_muts_unique = nt_muts.where(nt_muts == 0, 1).astype(int)
    total_nt_muts = int(nt_muts.values.sum())
    unique_nt_muts = int(nt_muts_unique.values.sum())

    nt_muts.reset_index(inplace=True)
    nt_muts_unique.reset_index(inplace=True)

    # Calculate mutation types
    trans_nt_muts = nt_muts.apply(
        lambda row: transversions_transitions(row['wt_nucleotide'], {nt: row[nt] for nt in 'ATGC'}),
        axis=1, result_type='expand'
    )
    trans_nt_muts_unique = nt_muts_unique.apply(
        lambda row: transversions_transitions(row['wt_nucleotide'], {nt: row[nt] for nt in 'ATGC'}),
        axis=1, result_type='expand'
    )
    all_mut_types = nt_muts.apply(
        lambda row: mut_type(row['wt_nucleotide'], {nt: row[nt] for nt in 'ATGC'}),
        axis=1, result_type='expand'
    )
    all_mut_types_unique = nt_muts_unique.apply(
        lambda row: mut_type(row['wt_nucleotide'], {nt: row[nt] for nt in 'ATGC'}),
        axis=1, result_type='expand'
    ).add_suffix('_unique')

    # Calculate basic NT statistics
    mean_nt_muts_per_seq = compute_mean_from_dist(nt_dist)
    total_insertion_length = df_dict['genotypes']['NT_insertion_length'].sum()
    total_deletion_length = df_dict['genotypes']['NT_deletion_length'].sum()

    values_list.extend([
        total_nt_muts, unique_nt_muts,
        mean_nt_muts_per_seq / reference_length,
        mean_nt_muts_per_seq,
        compute_median_from_dist(nt_dist),
        int(trans_nt_muts['transversions'].sum()),
        int(trans_nt_muts['transitions'].sum()),
        int(trans_nt_muts_unique['transversions'].sum()),
        int(trans_nt_muts_unique['transitions'].sum()),
        int(total_insertion_length),
        int(total_deletion_length)
    ])

    # Add individual mutation type counts (convert to int)
    values_list.extend([int(all_mut_types[mut_type].sum()) for mut_type in all_mut_types])
    values_list.extend([int(all_mut_types_unique[mut_type].sum()) for mut_type in all_mut_types_unique])

    return values_list


def build_column_names():
    """Build column names for statistics DataFrame."""
    cols = [
        'tag', 'barcode_group', 'total_seqs', 'total_unique_seqs', 'total_failed_seqs',
        'total_AA_mutations', 'unique_AA_mutations', 'mean_AA_mutations_per_seq',
        'median_AA_mutations_per_seq', 'total_NT_mutations', 'unique_NT_mutations',
        'mean_NT_mutations_per_base', 'mean_NT_mutations_per_seq',
        'median_NT_mutations_per_seq', 'total_transversions', 'total_transitions',
        'unique_transversions', 'unique_transitions', 'total_insertion_length',
        'total_deletion_length'
    ]

    # Add individual mutation type columns
    nts = 'ATGC'
    for wt in nts:
        for mut in nts:
            if wt != mut:
                cols.append(f'{wt}->{mut}')

    for wt in nts:
        for mut in nts:
            if wt != mut:
                cols.append(f'{wt}->{mut}_unique')

    return cols


def calculate_all_statistics(input_list, do_aa_analysis=False,
                             mutations_frequencies_raw=False):
    """
    Calculate statistics for all tag/barcode combinations.

    Args:
        input_list: List of input file paths
        do_aa_analysis: Whether AA mutation analysis was performed
        mutations_frequencies_raw: Whether frequencies are raw proportions (default: False)

    Returns:
        pandas.DataFrame: Statistics for all samples
    """
    stats_list = []
    file_dict = build_input_file_dict(input_list)

    for tag in file_dict:
        # Determine datatypes needed
        datatypes = ['genotypes', 'failures', 'NT-mutation-frequencies', 'NT-mutation-distribution']
        if do_aa_analysis:
            datatypes.extend(['AA-mutation-distribution', 'AA-mutation-frequencies'])

        for barcode_group in file_dict[tag]:
            # Load all data files
            df_dict = {dtype: pd.read_csv(file_dict[tag][barcode_group][dtype], index_col=0)
                      for dtype in datatypes}

            values = calculate_statistics_for_sample(
                df_dict, tag, barcode_group,
                do_aa_analysis, mutations_frequencies_raw
            )
            stats_list.append(values)

    # Build DataFrame
    cols = build_column_names()
    stats_df = pd.DataFrame(stats_list, columns=cols)
    stats_df.sort_values('barcode_group', inplace=True)

    # Round numeric columns
    stats_df['mean_NT_mutations_per_base'] = stats_df['mean_NT_mutations_per_base'].round(10)
    stats_df['mean_NT_mutations_per_seq'] = stats_df['mean_NT_mutations_per_seq'].round(2)

    return stats_df


def main():
    """Main entry point for mutation statistics calculation."""
    # Parse arguments
    if 'snakemake' in globals():
        # Running as snakemake script
        input_list = snakemake.input
        output_file = str(snakemake.output)
        do_aa_analysis = snakemake.params.do_aa_analysis
        mutations_frequencies_raw = snakemake.params.get('mutations_frequencies_raw', False)
    else:
        # Running as standalone script
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('--input', nargs='+', required=True,
                          help='Input mutation analysis files')
        parser.add_argument('--output', required=True,
                          help='Output statistics CSV file')
        parser.add_argument('--do-aa-analysis', action='store_true',
                          help='Whether AA mutation analysis was performed')
        parser.add_argument('--mutations-frequencies-raw', action='store_true',
                          help='Whether frequencies are raw proportions')
        args = parser.parse_args()

        input_list = args.input
        output_file = args.output
        do_aa_analysis = args.do_aa_analysis
        mutations_frequencies_raw = args.mutations_frequencies_raw

    # Calculate statistics
    stats_df = calculate_all_statistics(input_list, do_aa_analysis, mutations_frequencies_raw)

    # Save output
    stats_df.to_csv(output_file, index=False)


if __name__ == '__main__':
    main()
