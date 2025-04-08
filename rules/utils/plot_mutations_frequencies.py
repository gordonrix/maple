import numpy as np
import pandas as pd
import holoviews as hv
from common import colormaps, export_svg_plots, conspicuous_mutations
from snakemake.io import Namedlist
hv.extension("bokeh")

def main(frequencies_input, stats_input, output, labels, raw, num_positions, heatmap, colormap, export_SVG, NTorAA):
    mut_stats_df = pd.read_csv(stats_input, dtype={'tag':str,'barcode_group':str})

    if type(frequencies_input) != Namedlist:
        frequencies_input = [frequencies_input] # necessary for single input file

    if raw:
        value_label = 'total_count'
    else:
        value_label = 'proportion_of_seqs'

    if not heatmap:
        colormap = colormaps[NTorAA]

    plots = [] # for plotting all positions with numerical x axis
    plots_categorical_x = [] # to allow for slicing out most and least frequent positions
    ordered_columns = ['group', 'wt', 'position', 'mutation', 'total_count', 'proportion_of_seqs']
    dfs = []

    for i, in_file_full in enumerate(frequencies_input):

        # get mutation data, convert to tidy format
        in_file = in_file_full.split('/')[-1]
        tag = in_file.split('_')[-3]
        barcodeGroup = in_file.split('_')[-2]
        muts_df = pd.read_csv(in_file_full, index_col=False)
        muts_df = muts_df.rename(columns={muts_df.columns[0]:'wt'})
        muts_df = muts_df.melt(id_vars='wt', var_name='mutation', value_name=value_label)
        muts_df['position'] = pd.to_numeric(muts_df['wt'].str[1:])
        muts_df['wt'] = muts_df['wt'].str[0]
        plot_title = labels[i]
        muts_df['group'] = plot_title
        total_seqs = mut_stats_df.loc[(mut_stats_df['tag']==tag) & (mut_stats_df['barcode_group']==barcodeGroup), 'total_seqs'].iloc[0]

        if raw:
            muts_df['proportion_of_seqs'] = muts_df['total_count'] / total_seqs
        else:
            muts_df['total_count'] = muts_df['proportion_of_seqs'] * total_seqs

        dfs.append(muts_df)
        plots.append( conspicuous_mutations(muts_df, total_seqs, colormap=colormap, heatmap=heatmap).opts(title=plot_title, width=1000) )
        total_positions = len(muts_df['position'].unique())
        plots_categorical_x.append( conspicuous_mutations(muts_df, total_seqs, num_positions=total_positions, colormap=colormap, heatmap=heatmap).opts(title=plot_title, width=1000) )
    
    muts_df_all = pd.concat(dfs)[ordered_columns]

    # determine the most and least frequently mutated positions across all sample groups
    muts_df_all = muts_df_all.astype({'total_count':int, 'proportion_of_seqs':float})
    group_df = (muts_df_all.groupby(['wt', 'position'], as_index=False)
                            .sum()
                            .sort_values(['total_count'], ascending=True)
                            .reset_index()
                            .drop(columns='index'))
    group_df['wt_position'] = group_df['wt'] + group_df['position'].astype(str)
    least_frequent_positions = group_df.iloc[:num_positions]['wt_position'].to_list()
    most_frequent_positions = group_df.iloc[-num_positions:]['wt_position'].to_list()

    plots_least_frequent, plots_most_frequent = [],[]
    for plot in plots_categorical_x:
        plots_least_frequent.append(plot[least_frequent_positions,:].opts(width=num_positions*20))
        plots_most_frequent.append(plot[most_frequent_positions,:].opts(width=num_positions*20))

    hv.save( hv.Layout(plots_least_frequent).cols(1), output.least_frequent, backend='bokeh')
    hv.save( hv.Layout(plots_most_frequent).cols(1), output.most_frequent, backend='bokeh')
    hv.save( hv.Layout(plots).cols(1), output.all_muts, backend='bokeh')

    if export_SVG:
        export_svg_plots(plots, output.all_muts, labels=labels, export=export_SVG)

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
         snakemake.wildcards.NTorAA)

