import numpy as np
import pandas as pd
import holoviews as hv
import panel as pn
from common import colormaps, export_svg_plots, conspicuous_mutations
from snakemake.io import Namedlist
hv.extension("bokeh")

def get_top_references_by_abundance(mut_stats_df, tag, barcode_group, n_refs):
    """
    Get top N references by sequence abundance from mutation-stats DataFrame.

    Args:
        mut_stats_df: DataFrame from mutation-stats.csv with columns including
                      'tag', 'reference_name', 'barcode_group', 'total_seqs'
        tag: Tag string to filter by
        barcode_group: Barcode group string to filter by
        n_refs: Number of top references to return

    Returns:
        list: Reference names sorted by total_seqs (descending), limited to top n_refs
    """
    # Filter to specific tag and barcode_group
    filtered = mut_stats_df[
        (mut_stats_df['tag'] == tag) &
        (mut_stats_df['barcode_group'] == barcode_group)
    ].copy()

    # Sort by total_seqs and return top N reference names
    top_refs = filtered.nlargest(n_refs, 'total_seqs')['reference_name'].tolist()

    return top_refs

def main(frequencies_input, stats_input, output, labels, raw, num_positions, heatmap, colormap, export_SVG, NTorAA, split_by_reference):
    mut_stats_df = pd.read_csv(stats_input, dtype={'tag':str,'barcode_group':str})

    if type(frequencies_input) != Namedlist:
        frequencies_input = [frequencies_input] # necessary for single input file

    if raw:
        value_label = 'total_count'
    else:
        value_label = 'proportion_of_seqs'

    if not heatmap:
        colormap = colormaps[NTorAA]

    ordered_columns = ['group', 'wt', 'position', 'mutation', 'total_count', 'proportion_of_seqs']
    dfs = []

    # Collect all plots organized by [sample][reference]
    plots_by_sample_ref = {}  # {sample: {ref: (plot_numerical, plot_categorical)}}

    for i, in_file_full in enumerate(frequencies_input):

        in_file = in_file_full.split('/')[-1]
        tag = in_file.split('_')[-3]
        barcodeGroup = in_file.split('_')[-2]
        muts_df = pd.read_csv(in_file_full, index_col=False)
        # Keep wildtype entries (where wt == mutation, total_count is NaN) for plotting
        # The conspicuous_mutations function will handle them by setting proportion to -0.000001 for white coloring
        plot_title = labels[i]
        total_seqs = mut_stats_df.loc[(mut_stats_df['tag']==tag) & (mut_stats_df['barcode_group']==barcodeGroup), 'total_seqs'].iloc[0]

        # Get top N references by abundance from mutation-stats
        refs_to_plot = get_top_references_by_abundance(mut_stats_df, tag, barcodeGroup, split_by_reference)

        plots_by_sample_ref[plot_title] = {}

        # Create plots for each reference (will be arranged as columns)
        for ref in refs_to_plot:
            ref_df = muts_df[muts_df['reference_name'] == ref].copy()
            ref_df['group'] = plot_title
            dfs.append(ref_df)

            # Get total_seqs for this specific reference
            ref_total_seqs = mut_stats_df.loc[
                (mut_stats_df['tag'] == tag) &
                (mut_stats_df['barcode_group'] == barcodeGroup) &
                (mut_stats_df['reference_name'] == ref),
                'total_seqs'
            ].iloc[0]

            ref_title = f"partition: {barcodeGroup}\nreference: {ref}"

            plot_numerical = conspicuous_mutations(ref_df, int(ref_total_seqs), colormap=colormap, heatmap=heatmap, axis_type='numerical').opts(title=ref_title)
            total_positions = len(ref_df['position'].unique())
            plot_categorical = conspicuous_mutations(ref_df, int(ref_total_seqs), num_positions=total_positions, colormap=colormap, heatmap=heatmap).opts(title=ref_title, width=1000)

            plots_by_sample_ref[plot_title][ref] = (plot_numerical, plot_categorical)

    muts_df_all = pd.concat(dfs)[ordered_columns]

    # determine the most and least frequently mutated positions across all sample groups
    # Drop wildtype entries (NaN total_count) for aggregation purposes only
    muts_df_all = muts_df_all.dropna(subset=['total_count'])
    muts_df_all = muts_df_all.astype({'total_count':int, 'proportion_of_seqs':float})
    group_df = (muts_df_all.groupby(['wt', 'position'], as_index=False)
                            .sum()
                            .sort_values(['total_count'], ascending=True)
                            .reset_index()
                            .drop(columns='index'))
    group_df['wt_position'] = group_df['wt'] + group_df['position'].astype(str)
    least_frequent_positions = group_df.iloc[:num_positions]['wt_position'].to_list()
    most_frequent_positions = group_df.iloc[-num_positions:]['wt_position'].to_list()

    # Create Panel layouts: references as columns, samples as rows
    # Get all references (ordered by first sample's abundance)
    first_sample = list(plots_by_sample_ref.keys())[0]
    all_refs = list(plots_by_sample_ref[first_sample].keys())

    # Create columns for each reference
    columns_all = []
    columns_common = []
    columns_rare = []

    for ref in all_refs:
        # Collect plots for this reference across all samples
        ref_plots_all = []
        ref_plots_common = []
        ref_plots_rare = []

        for sample in plots_by_sample_ref.keys():
            if ref in plots_by_sample_ref[sample]:
                plot_numerical, plot_categorical = plots_by_sample_ref[sample][ref]
                ref_plots_all.append(plot_numerical)
                ref_plots_common.append(plot_categorical[most_frequent_positions,:].opts(width=num_positions*25))
                ref_plots_rare.append(plot_categorical[least_frequent_positions,:].opts(width=num_positions*25))

        # Create column for this reference
        columns_all.append(pn.Column(*ref_plots_all))
        columns_common.append(pn.Column(*ref_plots_common))
        columns_rare.append(pn.Column(*ref_plots_rare))

    # Arrange columns horizontally
    panel_layout_all = pn.Row(*columns_all)
    panel_layout_common = pn.Row(*columns_common)
    panel_layout_rare = pn.Row(*columns_rare)

    # Save using Panel
    panel_layout_rare.save(output.least_frequent)
    panel_layout_common.save(output.most_frequent)
    panel_layout_all.save(output.all_muts)

    if export_SVG:
        # Flatten plots for SVG export with corresponding labels
        all_plots = []
        svg_labels = []
        for sample in plots_by_sample_ref.keys():
            for ref in plots_by_sample_ref[sample].keys():
                plot_numerical, _ = plots_by_sample_ref[sample][ref]
                all_plots.append(plot_numerical)
                svg_labels.append(f"{sample}_{ref}")
        export_svg_plots(all_plots, output.all_muts, labels=svg_labels, export=export_SVG)

if __name__ == '__main__':
    main(snakemake.input.genotypes,
         snakemake.input.mut_stats,
         snakemake.output,
         snakemake.params.labels,
         snakemake.params.mutations_frequencies_raw,
         snakemake.params.number_of_positions,
         snakemake.params.heatmap,
         snakemake.params.colormap,
         snakemake.params.export_SVG,
         snakemake.wildcards.NTorAA,
         snakemake.params.split_by_reference)

