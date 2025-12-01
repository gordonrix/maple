#!/usr/bin/env python3
"""
plot_distribution.py

Script for maple pipeline. Generates bar plots and cumulative distribution plots.

Input data format:
    CSV files with the following required columns (order does not matter):
    - X-axis variable: What the distribution is counting, can be anything (e.g., 'mutations', 'hamming_distance')
    - Total counts: Must start with 'total ' (e.g., 'total sequences')
    - Proportion: Must start with 'proportion of ' (e.g., 'proportion of sequences')
    - Cumulative proportion: Must start with 'cumulative proportion of '
    - reference_name: Optional column for multi-reference data

    Column naming is strict - the script identifies columns by these name patterns.

    Example columns: ['mutations', 'total sequences', 'proportion of sequences',
                      'cumulative proportion of sequences', 'reference_name']

Output:
    - HTML file with interactive plots
    - Optional SVG exports for individual plots
    - Cumulative distribution overlay plot (all samples)
    - Individual distribution plots per sample/label
    - Grid layout for multi-reference data (rows=labels, columns=references)

Supports both snakemake mode and standalone command-line mode.

Author: Gordon Rix
"""

import pandas as pd
import numpy as np
import holoviews as hv
import hvplot.pandas
import panel as pn
import colorcet as cc
import argparse
import logging
from common import export_svg_plots, get_colors, str_to_bool

# Suppress bokeh validation warnings about fixed sizing mode
logging.getLogger('bokeh.core.validation.check').setLevel(logging.ERROR)

hv.extension('bokeh')

# Check if running in snakemake mode
try:
    snakemake
    snakemake_mode = True
except NameError:
    snakemake_mode = False

def main(input, output, labels, title, legendLabel, background, raw, export_svgs, cmap, x_range=False, y_range=False, split_by_reference=False, log_x=False, log_y=False, grid_rows_by_ref=False):
    """
    Main function that determines plotting layout based on input dimensions.

    Plotting modes:
        Case 1 (1D): Multiple labels, no reference splitting
                     → Rows = labels (first row cumulative overlay)
        Case 2 (1D): Single label, split by references
                     → Rows = references (first row cumulative overlay)
        Case 3 (2D): Multiple labels AND references
                     → Grid: rows = labels (first row cumulative overlay), columns = references (no cumulative)
    """
    x_range = convert_range(x_range)
    y_range = convert_range(y_range)
    inputDFs = [pd.read_csv(CSV, index_col=False) for CSV in input]

    # Add label column to each DataFrame and concatenate
    for df, label in zip(inputDFs, labels):
        df['label'] = label
    combined_df = pd.concat(inputDFs, ignore_index=True)

    # Determine plotting dimensions
    num_labels = len(labels)
    has_references = 'reference_name' in combined_df.columns
    # Check if any individual file has more than one reference
    has_multi_refs = any(len(df['reference_name'].unique()) > 1 for df in inputDFs) if has_references else False
    should_split_refs = has_references and split_by_reference is not False and has_multi_refs

    # Cases 1 and 2: 1D plotting (split by one dimension)
    if not should_split_refs:
        # Case 1: Split by labels only
        split_dim = 'label'
        split_values = labels
        dimension_label = legendLabel
    elif num_labels == 1:
        # Case 2: Split by references only
        split_dim = 'reference_name'
        split_values = get_filtered_references([combined_df], split_by_reference)
        dimension_label = 'reference'

    if not should_split_refs or num_labels == 1:
        # Handle Cases 1 and 2 with plot_1d
        colors = get_colors(split_values, cmap, background)
        plots, svg_labels = plot_1d(combined_df, split_dim, split_values, colors, title, legendLabel, raw, x_range, y_range, log_x, log_y, dimension_label)
        # Wrap in Panel Column for 1D layout
        panel_layout = pn.Column(*plots)
    else:
        # Case 3: Split by both labels AND references - 2D grid
        all_references = get_filtered_references([combined_df], split_by_reference)
        # Generate colors for whichever dimension will be in rows (for cumulative plot)
        row_values = all_references if grid_rows_by_ref else labels
        colors = get_colors(row_values, cmap, background)
        plots, svg_labels, panel_layout = plot_2d_grid(combined_df, labels, colors, all_references, title, legendLabel, raw, x_range, y_range, log_x, log_y, grid_rows_by_ref)

    # Export SVGs if requested
    if export_svgs:
        export_svg_plots(plots, output, labels=svg_labels, export=export_svgs, dimensions=(800, 600))

    # Save HTML output
    panel_layout.save(output)

def plot_1d(combined_df, split_dim, split_values, colors, title, legendLabel, raw, x_range, y_range, log_x, log_y, dimension_label=None, title_suffix=None):
    """
    Create 1D plot layout: cumulative + individual distribution plots in rows.

    Args:
        combined_df: Combined DataFrame with all data
        split_dim: Dimension column name to split by ('label' or 'reference_name')
        split_values: List of values to split by
        colors: List of colors for each split value
        dimension_label: Label for what dimension is being plotted (e.g., 'Sample', 'reference')

    Returns:
        tuple: (plots, svg_labels)
    """
    # Split the combined DataFrame by the dimension
    split_dfs = []
    for value in split_values:
        # Filter by split dimension
        df = combined_df[combined_df[split_dim] == value].copy()

        # Aggregate across other dimensions if needed
        if split_dim == 'label' and 'reference_name' in df.columns:
            # Case 1: Aggregating across references for each label
            df = aggregate_single_df(df)
        elif split_dim == 'reference_name':
            # Case 2: Drop label column
            df = df.drop(columns=['label'])

        split_dfs.append(df)

    # Create cumulative plot + individual plots
    plots = [plot_cumsum(split_dfs, split_values, colors, title, legendLabel, x_range=x_range, y_range=y_range, log_x=log_x, log_y=log_y, title_suffix=title_suffix)]

    # Use dimension_label if provided, otherwise use legendLabel
    if dimension_label is None:
        dimension_label = legendLabel

    # Build plot titles with optional suffix
    for df, label, c in zip(split_dfs, split_values, colors):
        plot_title = f"{dimension_label}: {label}"
        if title_suffix:
            plot_title += f"\n{title_suffix}"
        plots.append(plot_dist(df, color=c, title=plot_title, raw=raw, x_range=x_range, y_range=y_range, log_x=log_x, log_y=log_y))

    svg_labels = ['cumulative'] + list(split_values)
    return plots, svg_labels

def plot_2d_grid(combined_df, labels, colors, all_references, title, legendLabel, raw, x_range, y_range, log_x, log_y, grid_rows_by_ref=False):
    """
    Create 2D grid.

    Default (grid_rows_by_ref=False): rows = (cumulative + samples), columns = references
    If grid_rows_by_ref=True: rows = (cumulative + references), columns = samples

    Uses plot_1d for each column, then combines them horizontally.

    Returns:
        tuple: (flat_plots, svg_labels, panel_layout)
    """
    columns = []
    all_plots = []
    all_svg_labels = []

    # Determine which dimension to iterate over (columns) and which to split by (rows)
    if grid_rows_by_ref:
        outer_values = labels
        outer_dim = 'label'
        inner_values = all_references
        inner_dim = 'reference_name'
        outer_dim_label = legendLabel
        inner_dim_label = 'reference'
    else:
        outer_values = all_references
        outer_dim = 'reference_name'
        inner_values = labels
        inner_dim = 'label'
        outer_dim_label = 'reference'
        inner_dim_label = legendLabel

    # Create one column per outer dimension value
    for i, outer_val in enumerate(outer_values):
        # Filter combined_df to this outer dimension value
        filtered_df = combined_df[combined_df[outer_dim] == outer_val].copy()

        # Colors are always for the inner dimension (rows) - use same colors for all columns
        col_colors = colors

        # Call plot_1d for this column
        col_plots, col_svg_labels = plot_1d(filtered_df, inner_dim, inner_values, col_colors, title, legendLabel, raw, x_range, y_range, log_x, log_y, dimension_label=inner_dim_label, title_suffix=f"{outer_dim_label}: {outer_val}")

        # Add outer dimension suffix to svg labels
        col_svg_labels = [f"{label}_{outer_val}" for label in col_svg_labels]

        # Create Panel Column
        columns.append(pn.Column(*col_plots))
        all_plots.extend(col_plots)
        all_svg_labels.extend(col_svg_labels)

    # Combine columns horizontally
    panel_layout = pn.Row(*columns)

    return all_plots, all_svg_labels, panel_layout

def aggregate_single_df(df):
    """
    Aggregate a single DataFrame across references by summing.
    Recalculates proportions and cumulative sums after aggregation.

    Args:
        df: DataFrame with reference_name column

    Returns:
        DataFrame with reference_name aggregated away
    """
    x_col = get_x_column(df)
    total_col = get_total_column(df)
    proportion_col = get_proportion_column(df)
    cumsum_col = get_cumsum_column(df)

    # Group by x and sum the total column
    agg_df = df.groupby(x_col, as_index=False)[total_col].sum()

    # Recalculate proportions and cumulative sums
    total_sum = agg_df[total_col].sum()
    agg_df[proportion_col] = agg_df[total_col] / total_sum
    agg_df[cumsum_col] = agg_df[proportion_col].cumsum()

    return agg_df

def get_filtered_references(inputDFs, split_by_reference):
    """Get list of references to plot, filtered by abundance if split_by_reference is an int."""
    all_references = sorted(set([ref for df in inputDFs for ref in df['reference_name'].unique()]))

    # Check for int but not bool (since bool is subclass of int in Python)
    if isinstance(split_by_reference, int) and not isinstance(split_by_reference, bool):
        # Calculate total abundance per reference across all files
        ref_abundance = {}
        for df in inputDFs:
            total_col = get_total_column(df)
            for ref in df['reference_name'].unique():
                ref_df = df[df['reference_name'] == ref]
                total = ref_df[total_col].sum()
                ref_abundance[ref] = ref_abundance.get(ref, 0) + total
        # Sort and take top N
        all_references = sorted(ref_abundance.keys(), key=lambda x: ref_abundance[x], reverse=True)[:split_by_reference]

    return all_references

def get_x_column(df):
    """Get the x-axis column name."""
    return [col for col in df.columns if col not in ['reference_name', 'label'] and
            not col.startswith('total') and not col.startswith('proportion') and
            not col.startswith('cumulative')][0]

def get_total_column(df):
    """Get the total counts column name."""
    return [col for col in df.columns if col.startswith('total ')][0]

def get_proportion_column(df):
    """Get the proportion column name."""
    return [col for col in df.columns if col.startswith('proportion of ')][0]

def get_cumsum_column(df):
    """Get the cumulative proportion column name."""
    return [col for col in df.columns if col.startswith('cumulative proportion of ')][0]

def convert_range(range):
    """
    converts a range string to a tuple of floats, or None for blank (e.g. ',10' becomes (None, 100))
    """
    if range:
        range_list = range.split(',')
        range_list = [float(val) if val else None for val in range_list]
        return tuple(range_list)
    else:
        return (None, None)

def plot_cumsum(DFlist, labels, colors, title, legendLabel, x_range=False, y_range=False, log_x=False, log_y=False, title_suffix=None):
    """
    generates a distribution hv.Curve plot showing cumulative sums for all input files as hv.Curve on an hv.Overlay plot

    DFlist:         list of pd.DataFrames that contain cumulative sums in the format generated by the common.dist_to_DF function
    labels:         list of strings, names to use to label each group for cumulative distribution hv.Overlay
    colors:         list of strings, color values index matched to DFlist to be used to color each line
    title:          string, label for all samples to include in the title
    legendLabel:    string, label to add to the legend (only appears in hover)
    x_max:          int, maximum value for x axis
    y_max:          int, maximum value for y axis
    """

    plotList = []
    for dfOG, label, color in zip(DFlist, labels, colors):

        df = dfOG.copy()    # prevent adding a column to the global dataframe

        # Find columns by name pattern
        x = get_x_column(df)
        total = get_total_column(df)
        proportion = get_proportion_column(df)
        cumsum = get_cumsum_column(df)

        y = total.split('total ')[1]
        df[legendLabel] = label

        plotList.append( hv.Curve(df, kdims=[x],
                            vdims = [cumsum, proportion, total] + [legendLabel], label=label).opts(
                                        hv.opts.Curve(line_width=3, tools=['hover'], color=color)) )

    # Calculate legend columns based on number of items (16 items per column)
    num_items = len(plotList)
    legend_cols = (num_items + 15) // 16  # 1 col for ≤16 items, 2 cols for 17-32, etc.

    # Build title with optional suffix
    plot_title = f"{x} among {y}, cumulative distribution\nsample: {title}"
    if title_suffix:
        plot_title += f"\n{title_suffix}"

    plot = hv.Overlay(plotList).opts(show_legend=True, legend_position='right', legend_cols=legend_cols,
                            frame_width=550, frame_height=500,
                            title=plot_title,
                            # legend_title = legendLabel, # This doesn't work for some reason, but a legend title would be nice
                            xlabel=x, ylabel='cumulative frequency', ylim=y_range, xlim=x_range,
                            logx=log_x, logy=log_y,
                            fontsize={'title':18,'labels':18,'xticks':16,'yticks':16, 'legend': 16})

    return plot


def plot_dist(df, color='grey', title='', raw=False, x_range=None, y_range=None, log_x=False, log_y=False):
    """
    generates a histogram with bin size of 1 for a distribution input file

    df:     pd.DataFrame, distribution data in the format generated by the common.dist_to_DF function
    color:  string, color to use for bars, all will be the same color
    title:  string, title that describes the individual sample
    raw:    bool, determines if y axis will show raw counts or proportion of the total
                default behavior is to show proportion of the total, not raw counts
    """

    # Find columns by name pattern
    x = get_x_column(df)
    total = get_total_column(df)
    proportion = get_proportion_column(df)
    cumsum = get_cumsum_column(df)

    y = total.split('total ')[1]

    totalSum = int(df[total].sum())

    if raw:
        vdims = [total, proportion, cumsum]
        ylabel = f"frequency as total"
    else:
        vdims = [proportion, cumsum, total]
        ylabel = f"frequency as proportion"

    ylabel += f'\n(n = {totalSum} {y})'

    hvDF = hv.Dataset(df, kdims=[x], vdims=vdims)
    plot = hv.Histogram(hvDF).opts(hv.opts.Histogram(tools=['hover'], color=color, frame_width=550, frame_height=500,
                                    title=f"{x} among {y}\n{title}",
                                    ylabel=ylabel, ylim=y_range, xlim=x_range,
                                    logx=log_x, logy=log_y,
                                    fontsize={'title':18,'labels':18,'xticks':16,'yticks':16, 'legend': 16}))

    return plot

if __name__ == '__main__':
    if snakemake_mode:
        log_x = getattr(snakemake.params, 'log_x', False)
        log_y = getattr(snakemake.params, 'log_y', False)
        grid_rows_by_ref = getattr(snakemake.params, 'grid_rows_by_ref', False)
        main(snakemake.input, snakemake.output.plot, snakemake.params.labels, snakemake.params.title,
             snakemake.params.legend_label, snakemake.params.background, snakemake.params.raw,
             snakemake.params.export_SVG, snakemake.params.colormap, snakemake.params.x_max,
             snakemake.params.y_max, snakemake.params.split_by_reference, log_x, log_y, grid_rows_by_ref)
    else:
        parser = argparse.ArgumentParser(description="Generate distribution plots from CSV files")
        parser.add_argument("--input", nargs='+', required=True, help="Input CSV file(s)")
        parser.add_argument("--output", required=True, help="Output HTML file path")
        parser.add_argument("--labels", nargs='+', required=True, help="Labels for each input file")
        parser.add_argument("--title", required=True, help="Title for the plots")
        parser.add_argument("--legend_label", required=True, help="Label for legend")
        parser.add_argument("--background", default='False', help="Background sample label or 'False'")
        parser.add_argument("--raw", action='store_true', help="Use raw counts instead of proportions")
        parser.add_argument("--export_svgs", default='False', help="Export SVGs: 'True', 'False', or filter string")
        parser.add_argument("--colormap", default='kbc_r', help="Colormap name")
        parser.add_argument("--x_range", default='False', help="X-axis range as 'min,max' or 'False'")
        parser.add_argument("--y_range", default='False', help="Y-axis range as 'min,max' or 'False'")
        parser.add_argument("--split_by_reference", default='False', help="Split by reference: 'False', 'True', or int N")
        parser.add_argument("--log_x", action='store_true', help="Use logarithmic scale for x-axis")
        parser.add_argument("--log_y", action='store_true', help="Use logarithmic scale for y-axis")
        parser.add_argument("--grid_rows_by_ref", action='store_true', help="Transpose 2D grid: rows by reference, columns by sample")

        args = parser.parse_args()

        # Convert background - can be False or a string label
        try:
            background = False if str_to_bool(args.background) == False else args.background
        except ValueError:
            background = args.background  # It's a label string, not a boolean

        # raw is already a boolean from action='store_true'
        raw = args.raw

        # Convert export_svgs - can be bool or string filter
        try:
            export_svgs = str_to_bool(args.export_svgs)
        except ValueError:
            export_svgs = args.export_svgs  # It's a filter string, not a boolean

        # Convert split_by_reference - can be False, True, or int
        try:
            split_by_reference = str_to_bool(args.split_by_reference)
        except ValueError:
            try:
                split_by_reference = int(args.split_by_reference)
            except ValueError:
                parser.error(f"--split_by_reference must be 'True', 'False', or an integer, got: {args.split_by_reference}")

        # Convert ranges - can be False or string like "0,100"
        try:
            x_range = False if str_to_bool(args.x_range) == False else args.x_range
        except ValueError:
            x_range = args.x_range  # It's a range string, not a boolean

        try:
            y_range = False if str_to_bool(args.y_range) == False else args.y_range
        except ValueError:
            y_range = args.y_range  # It's a range string, not a boolean

        main(args.input, args.output, args.labels, args.title, args.legend_label,
             background, raw, export_svgs, args.colormap, x_range, y_range, split_by_reference, args.log_x, args.log_y, args.grid_rows_by_ref)