#
#  DESCRIPTION   : Script for maple pipeline. Generates bar plots that describe the demultiplexing outcomes,
#                   including plotting for label_only barcodes. Supports multi-reference stacked bars.
#
#  AUTHOR(S)     : Gordon Rix
#

import argparse
import json
import pandas as pd
import numpy as np
import colorcet as cc
import holoviews as hv
from natsort import index_natsorted
from common import dist_to_DF, get_colors, str_to_bool
from plot_distribution import plot_cumsum, plot_dist
import hvplot.pandas

# Check if running in snakemake mode
try:
    snakemake
    snakemake_mode = True
except NameError:
    snakemake_mode = False


def get_filtered_references(df, split_by_reference):
    """
    Get list of top N references by abundance, where N = split_by_reference.
    Adapted from plot_distribution.py

    Args:
        df: DataFrame with 'reference_name' and 'demuxed_count' columns
        split_by_reference: int N to get top N references, or True for all

    Returns:
        list of reference names (top N by abundance)
    """
    # Calculate total counts per reference
    ref_abundance = df.groupby('reference_name')['demuxed_count'].sum().to_dict()

    # Sort by abundance and get all references
    all_references = sorted(ref_abundance.keys(), key=lambda x: ref_abundance[x], reverse=True)

    # Check for int but not bool (since bool is subclass of int in Python)
    if isinstance(split_by_reference, int) and not isinstance(split_by_reference, bool):
        # Return top N references
        return all_references[:split_by_reference]
    else:
        # Return all references if True
        return all_references


def sort_by_barcode_groups(df, bc_groups):
    """
    Sort dataframe by output_file_barcodes: user-provided names first, then natural sort.

    Args:
        df: DataFrame with 'output_file_barcodes' column
        bc_groups: List of user-provided barcode group names (order matters)

    Returns:
        Sorted DataFrame
    """
    df_names = df[df['output_file_barcodes'].isin(bc_groups)]
    df_names = df_names.sort_values(by='output_file_barcodes', key=lambda col: [bc_groups.index(x) for x in col])
    df_no_names = df[~df['output_file_barcodes'].isin(bc_groups)]
    df_no_names = df_no_names.sort_values(by='output_file_barcodes', key=lambda x: np.argsort(index_natsorted(df_no_names['output_file_barcodes'])))
    return pd.concat([df_names, df_no_names])

def main(input, output, barcode_info_dict, barcode_groups_dict, split_by_reference=False):
    # Load original data - keep this intact (never mutate) for label_only analysis
    demux_df_original = pd.read_csv(input)
    assert len(demux_df_original['tag'].unique())==1   # this script uses demux stats output for only a single tag
    demux_df_original = demux_df_original.drop(columns='tag')
    demux_df_original['unique_label_only_count'] = 1 # counter for unique label_only barcodes, if they are present

    color = 'grey'

    label_only_bc_types = [bc_type for bc_type in barcode_info_dict.keys() if barcode_info_dict[bc_type].get('label_only', False)]
    group_bc_types = [bc_type for bc_type in barcode_info_dict.keys() if bc_type not in label_only_bc_types]

    bc_groups = list(barcode_groups_dict.keys()) if barcode_groups_dict else []

    # Check if we have reference data and should split by reference
    has_references = 'reference_name' in demux_df_original.columns
    has_multi_refs = len(demux_df_original['reference_name'].unique()) > 1 if has_references else False
    should_split_refs = has_references and split_by_reference is not False and has_multi_refs

    # ========== PREPARE DATA FOR MAIN PLOT ==========
    # If splitting by reference, add 'reference_category' to original data (both plots need it)
    if should_split_refs:
        top_references = get_filtered_references(demux_df_original, split_by_reference)
        demux_df_original['reference_category'] = demux_df_original['reference_name'].apply(
            lambda x: x if x in top_references else 'other'
        )

    # Group data for main plot
    # Note: barcodes_count is the actual count for each (barcode_group, reference, label_only_variation)
    # demuxed_count is the total for the barcode group (repeated across rows)
    if should_split_refs:
        # Group by output_file_barcodes AND reference_category for stacked bars
        # Sum barcodes_count to get the actual count per reference
        aggs = {'barcodes_count': 'sum', 'unique_label_only_count': 'sum'}
        aggs.update({bc_type: 'first' for bc_type in group_bc_types})
        df_grouped = demux_df_original.groupby(['output_file_barcodes', 'reference_category'], as_index=False).agg(aggs)
        df_grouped = df_grouped.rename(columns={'barcodes_count': 'total_barcodes_count'})
    elif has_references:
        # Aggregate across references - first collapse label_only, then sum across references
        # Step 1: Group by (output_file_barcodes, reference_name) to collapse label_only variations
        temp_aggs = {'demuxed_count': 'first', 'unique_label_only_count': 'sum'}
        temp_aggs.update({bc_type: 'first' for bc_type in group_bc_types})
        df_temp = demux_df_original.groupby(['output_file_barcodes', 'reference_name'], as_index=False).agg(temp_aggs)
        # Step 2: Sum across references
        aggs = {'demuxed_count': 'sum', 'unique_label_only_count': 'sum'}
        aggs.update({bc_type: 'first' for bc_type in group_bc_types})
        df_grouped = df_temp.groupby('output_file_barcodes', as_index=False).agg(aggs)
        df_grouped = df_grouped.rename(columns={'demuxed_count': 'total_barcodes_count'})
    else:
        # No references - simple grouping, use 'first' since value is repeated per label_only variation
        aggs = {'output_file_barcodes': 'first', 'demuxed_count': 'first', 'unique_label_only_count': 'sum'}
        aggs.update({bc_type: 'first' for bc_type in group_bc_types})
        df_grouped = demux_df_original.groupby('output_file_barcodes', as_index=False).agg(aggs)
        df_grouped = df_grouped.rename(columns={'demuxed_count': 'total_barcodes_count'})

    # Sort by barcode groups (user-provided names first, then natural sort)
    df_grouped = sort_by_barcode_groups(df_grouped, bc_groups)

    # Create main demux plot (stacked by reference if applicable)
    if should_split_refs:
        # Get colors for reference categories (top N + 'other')
        ref_categories = top_references + ['other']
        ref_colors = get_colors(ref_categories, cmap='kbc_r', background=False)
        plot_color = ref_colors
        plot_by = 'reference_category'
        plot_stacked = True
    else:
        plot_color = color
        plot_by = None
        plot_stacked = False

    # Calculate total count for y-axis label
    total_count = int(df_grouped['total_barcodes_count'].sum())

    plot = df_grouped.hvplot.bar(
        x='output_file_barcodes',
        y='total_barcodes_count',
        by=plot_by,
        stacked=plot_stacked,
        hover_cols=[col for col in ['fwd', 'rvs'] if col in df_grouped.columns],
        color=plot_color,
        title="total demultiplex counts for each group",
        ylabel=f'total barcodes count (n = {total_count} sequences)',
        rot=70,
        fontsize={'title': 16, 'labels': 14, 'xticks': 10, 'yticks': 10},
        width=800,
        height=600,
        legend='top_right'
    )

    if label_only_bc_types:
        # ========== PREPARE DATA FOR LABEL_ONLY PLOTS ==========
        # Label_only bar plots use same reference splitting as main plot (stacked if enabled)
        # Only histograms aggregate across references (because they show distributions)

        # find all rows in which any label_only barcodes failed identification
        query = " | ".join([x + " == 'fail'" for x in label_only_bc_types])

        # Group label_only data (with reference splitting if enabled)
        # Note: barcodes_count is the actual count for each (barcode_group, reference, label_only_variation)
        if should_split_refs:
            # Group by output_file_barcodes AND reference_category for stacked bars
            # Get failed counts - sum barcodes_count for rows where label_only barcodes failed
            df_failed_grouped = demux_df_original.query(query).groupby(['output_file_barcodes', 'reference_category'], as_index=False).agg({
                'barcodes_count': 'sum'
            }).rename(columns={'barcodes_count': 'failed_label_only_barcodes_count'})

            # Get total counts - sum barcodes_count to get actual count per reference
            label_only_df = demux_df_original.groupby(['output_file_barcodes', 'reference_category'], as_index=False).agg({
                'barcodes_count': 'sum',
                'unique_label_only_count': 'sum'
            }).rename(columns={'barcodes_count': 'total_barcodes_count'})
        elif has_references:
            # Aggregate across references
            # Get failed counts - sum barcodes_count for failed label_only barcodes, then sum across refs
            df_failed_grouped = demux_df_original.query(query).groupby('output_file_barcodes', as_index=False).agg({
                'barcodes_count': 'sum'
            }).rename(columns={'barcodes_count': 'failed_label_only_barcodes_count'})

            # Total counts - first collapse label_only, then sum across references
            df_temp = demux_df_original.groupby(['output_file_barcodes', 'reference_name'], as_index=False).agg({
                'demuxed_count': 'first',
                'unique_label_only_count': 'sum'
            })
            label_only_df = df_temp.groupby('output_file_barcodes', as_index=False).agg({
                'demuxed_count': 'sum',
                'unique_label_only_count': 'sum'
            }).rename(columns={'demuxed_count': 'total_barcodes_count'})
        else:
            # No references
            df_failed_grouped = demux_df_original.query(query).groupby('output_file_barcodes', as_index=False).agg({
                'barcodes_count': 'sum'
            }).rename(columns={'barcodes_count': 'failed_label_only_barcodes_count'})

            label_only_df = demux_df_original.groupby('output_file_barcodes', as_index=False).agg({
                'demuxed_count': 'first',
                'unique_label_only_count': 'sum'
            }).rename(columns={'demuxed_count': 'total_barcodes_count'})

        # Merge failed counts into label_only_df
        merge_on = ['output_file_barcodes', 'reference_category'] if should_split_refs else ['output_file_barcodes']
        label_only_df = label_only_df.merge(df_failed_grouped, how='left', on=merge_on)
        label_only_df['failed_label_only_barcodes_count'] = label_only_df['failed_label_only_barcodes_count'].fillna(0).astype(int)
        label_only_df['success_label_only_barcodes_count'] = label_only_df['total_barcodes_count'] - label_only_df['failed_label_only_barcodes_count']

        # Sort label_only_df same way as df_grouped
        label_only_df = sort_by_barcode_groups(label_only_df, bc_groups)

        # plot number of sequences in each group with one or more failed label_only barcode identification (only if there are failures)
        if label_only_df['failed_label_only_barcodes_count'].sum() > 0:
            total_failed = int(label_only_df['failed_label_only_barcodes_count'].sum())
            plot = plot + label_only_df.hvplot.bar(
                x='output_file_barcodes',
                y='failed_label_only_barcodes_count',
                by=plot_by,
                stacked=plot_stacked,
                hover_cols=[col for col in ['fwd', 'rvs'] if col in label_only_df.columns],
                color=plot_color,
                title="failed label_only demultiplex counts for each group",
                ylabel=f'failed label only barcodes count (n = {total_failed} sequences)',
                rot=70,
                fontsize={'title': 16, 'labels': 14, 'xticks': 12, 'yticks': 12},
                width=800,
                height=600,
                legend='top_right'
            )

        # plot number of sequences in each group with successful identification of all label_only barcodes
        if label_only_df['success_label_only_barcodes_count'].sum() > 0:
            total_success = int(label_only_df['success_label_only_barcodes_count'].sum())
            plot = plot + label_only_df.hvplot.bar(
                x='output_file_barcodes',
                y='success_label_only_barcodes_count',
                by=plot_by,
                stacked=plot_stacked,
                hover_cols=[col for col in ['fwd', 'rvs'] if col in label_only_df.columns],
                color=plot_color,
                title="successful label_only demultiplex counts for each group",
                ylabel=f'success label only barcodes count (n = {total_success} sequences)',
                rot=70,
                fontsize={'title': 16, 'labels': 14, 'xticks': 12, 'yticks': 12},
                width=800,
                height=600,
                legend='top_right'
            )

        # plot number of unique label_only barcodes in each group
        if label_only_df['unique_label_only_count'].sum() > 0:
            total_unique = int(label_only_df['unique_label_only_count'].sum())
            plot = plot + label_only_df.hvplot.bar(
                x='output_file_barcodes',
                y='unique_label_only_count',
                by=plot_by,
                stacked=plot_stacked,
                hover_cols=[col for col in ['fwd', 'rvs'] if col in label_only_df.columns],
                color=plot_color,
                title="unique label_only demultiplex counts for each group",
                ylabel=f'unique label only count (n = {total_unique} unique barcodes)',
                rot=70,
                fontsize={'title': 16, 'labels': 14, 'xticks': 12, 'yticks': 12},
                width=800,
                height=600,
                legend='top_right'
            )

        # get distributions of label_only barcode counts for each barcode group, then plot
        # Histograms aggregate across references (showing distribution of barcodes_count)
        hist_list = []
        no_fail = demux_df_original.loc[~demux_df_original.eval(query)]

        # Loop through output_file_barcodes in same sorted order as label_only_df
        for output_file_barcode in label_only_df['output_file_barcodes'].unique():
            if bc_groups:
                if output_file_barcode not in bc_groups:
                    continue
            subset = no_fail[no_fail['output_file_barcodes'] == output_file_barcode]
            if len(subset) == 0:
                continue
            max_x_val = int(subset['barcodes_count'].max())
            total_unique_groups = len(subset)

            hist_list.append(subset.hvplot.hist(
                'barcodes_count',
                bins=np.arange(0.5, max_x_val + 1.5, 1),
                title=f"{output_file_barcode}: label_only demux count distribution",
                color='grey',
                xlabel='label_only demux count',
                ylabel=f'unique barcode groups (n = {total_unique_groups} groups)',
                fontsize={'title': 16, 'labels': 14, 'xticks': 12, 'yticks': 12},
                width=800,
                height=600
            ))

        if hist_list:
            hist_plots = hv.Layout(hist_list)
            plot = plot + hist_plots
            plot = plot.cols(1).opts(shared_axes=False)

    # Reference coverage distribution (only if we have multiple references)
    if has_multi_refs:
        ref_coverage_list = []

        # Loop through barcode groups to create one plot per group
        for output_file_barcode in df_grouped['output_file_barcodes'].unique() if not should_split_refs else label_only_df['output_file_barcodes'].unique():
            if bc_groups:
                if output_file_barcode not in bc_groups:
                    continue

            # Get reference counts for this barcode group (sum across label_only variations)
            subset = demux_df_original[demux_df_original['output_file_barcodes'] == output_file_barcode]
            ref_counts = subset.groupby('reference_name')['barcodes_count'].sum()

            if len(ref_counts) == 0:
                continue

            # Calculate log10(count) bins: 0 for count=1, 1 for count=2-10, 2 for count=11-100, etc.
            log_bins = []
            for count in ref_counts:
                if count == 1:
                    log_bins.append(0)
                else:
                    log_bins.append(int(np.floor(np.log10(count))))

            log_df = pd.DataFrame({'log_bin': log_bins})
            max_log_bin = log_df['log_bin'].max()
            total_refs = len(ref_counts)

            # Create histogram with integer bins
            ref_coverage_list.append(log_df.hvplot.hist(
                'log_bin',
                bins=np.arange(-0.5, max_log_bin + 1.5, 1),
                title=f"{output_file_barcode}: reference coverage distribution",
                color='grey',
                xlabel='reference coverage (log10(max of bin))',
                ylabel=f'number of references (n = {total_refs} references)',
                fontsize={'title': 16, 'labels': 14, 'xticks': 12, 'yticks': 12},
                width=800,
                height=600
            ))

        if ref_coverage_list:
            ref_coverage_plots = hv.Layout(ref_coverage_list)
            plot = plot + ref_coverage_plots
            plot = plot.cols(1).opts(shared_axes=False)

    hvplot.save(plot, output)

if __name__ == '__main__':
    if snakemake_mode:
        split_by_reference = getattr(snakemake.params, 'distribution_split_by_reference', False)
        main(snakemake.input.CSV, snakemake.output.plot, snakemake.params.barcodeInfo,
             snakemake.params.barcodeGroups, split_by_reference)
    else:
        parser = argparse.ArgumentParser(description="Generate demultiplex bar plots from CSV file")
        parser.add_argument("--input", required=True, help="Input demux stats CSV file")
        parser.add_argument("--output", required=True, help="Output HTML file path")
        parser.add_argument("--barcode_info", required=True, help="Barcode info dict (as JSON string)")
        parser.add_argument("--barcode_groups", required=True, help="Barcode groups dict (as JSON string)")
        parser.add_argument("--split_by_reference", default='False',
                          help="Split by reference: 'False', 'True', or int N for top N references")

        args = parser.parse_args()

        # Convert barcode_info and barcode_groups from JSON strings
        barcode_info_dict = json.loads(args.barcode_info)
        barcode_groups_dict = json.loads(args.barcode_groups)

        # Convert split_by_reference - can be False, True, or int
        try:
            split_by_reference = str_to_bool(args.split_by_reference)
        except ValueError:
            try:
                split_by_reference = int(args.split_by_reference)
            except ValueError:
                parser.error(f"--split_by_reference must be 'True', 'False', or an integer, got: {args.split_by_reference}")

        main(args.input, args.output, barcode_info_dict, barcode_groups_dict, split_by_reference)