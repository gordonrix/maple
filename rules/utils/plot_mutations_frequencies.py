import numpy as np
import pandas as pd
import holoviews as hv
from common import colormaps, export_svg_plots, conspicuous_mutations
from snakemake.io import Namedlist
hv.extension("bokeh")


def main(frequencies_input, stats_input, output, mutations_frequencies_raw, heatmap, colormap, export_SVG, wildcards):
    number_of_positions = 20    # number of positions to include for most and least frequent

    mutStatsDF = pd.read_csv(stats_input, dtype={'tag':str,'barcode_group':str})

    if type(frequencies_input) != Namedlist:
        frequencies_input = [frequencies_input]

    if mutations_frequencies_raw:
        value_label = 'total_count'
    else:
        value_label = 'proportion_of_seqs'

    if not heatmap:
        colormap = colormaps[wildcards.NTorAA]

    plots = []
    mutsDF_all = pd.DataFrame(columns=['group', 'wt', 'position', 'mutation', 'total_count'])

    group_list = []
    for inFile in frequencies_input:

        # Get mutation data, convert to tidy format
        barcodeGroup = inFile.split('_')[-2]
        group_list.append(barcodeGroup)
        plotTitle = f'{wildcards.tag}_{barcodeGroup}'
        mutsDF = pd.read_csv(inFile, index_col=False)
        mutsDF = mutsDF.rename(columns={mutsDF.columns[0]:'wt'})
        wt = list(mutsDF.columns[1:])
        total_seqs = mutStatsDF.loc[mutStatsDF['barcode_group']==barcodeGroup, 'total_seqs'].iloc[0]
        mutsDF = mutsDF.melt(id_vars='wt', var_name='mutation', value_name=value_label)
        mutsDF['position'] = pd.to_numeric(mutsDF['wt'].str[1:])
        mutsDF['wt'] = mutsDF['wt'].str[0]
        mutsDF['group']=plotTitle

        if mutations_frequencies_raw:
            mutsDF['proportion_of_seqs'] = mutsDF['total_count'] / total_seqs
        else:
            mutsDF['total_count'] = mutsDF['proportion_of_seqs'] * total_seqs

        plots.append( conspicuous_mutations(mutsDF, total_seqs, colormap=colormap, heatmap=heatmap).opts(title=plotTitle, width=1000) )
        mutsDF_all = pd.concat([mutsDF_all, mutsDF])

    # determine the most and least frequently mutated positions across all sample groups
    mutsDF_all = mutsDF_all.astype({'total_count':int, 'proportion_of_seqs':float})
    group_DF = mutsDF_all.groupby(['wt', 'position'], as_index=False).sum().sort_values(['total_count'], ascending=True).reset_index().drop(columns='index')
    group_DF['wt_position'] = group_DF['wt'] + group_DF['position'].astype(str)
    least_frequent_positions = group_DF.iloc[:number_of_positions]['wt_position'].to_list()
    most_frequent_positions = group_DF.iloc[-number_of_positions:]['wt_position'].to_list()

    plots_least_frequent, plots_most_frequent = [],[]
    for plot in plots:
        plots_least_frequent.append(plot[least_frequent_positions,:].opts(width=number_of_positions*20))
        plots_most_frequent.append(plot[most_frequent_positions,:].opts(width=number_of_positions*20))

    hv.save( hv.Layout(plots_least_frequent).cols(1), output.least_frequent, backend='bokeh')
    hv.save( hv.Layout(plots_most_frequent).cols(1), output.most_frequent, backend='bokeh')
    hv.save( hv.Layout(plots).cols(1), output.all_muts, backend='bokeh')

    if export_SVG:
        export_svg_plots(plots, output.all_muts, labels=group_list, export=export_SVG)

if __name__ == '__main__':
    main(snakemake.input.genotypes, snakemake.input.mut_stats, snakemake.output, snakemake.params.mutations_frequencies_raw, snakemake.params.heatmap, snakemake.params.cmap, snakemake.params.export_SVG, snakemake.wildcards)

