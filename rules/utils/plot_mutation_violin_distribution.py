#
#  DESCRIPTION   : Script for maple pipeline. Generates violin plots for NT and (optionally) AA mutations.
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import holoviews as hv
import hvplot.pandas
import colorcet as cc
from itertools import product
from common import export_svg_plots, get_colors, cmap_dict
hv.extension('bokeh')

def main(input, output, group_col, x_label, export_svgs, cmap, background):

    df = pd.read_csv(input, index_col=False)
    include_AA = True if 'AA_substitutions_nonsynonymous_count' in df.columns else False

    plots = plot_violin(df, include_AA, group_col, x_label, cmap, background)

    if export_svgs:
        export_svg_plots(plots, output, export=export_svgs)
    plots = hv.Layout(plots).cols(1)
    hvplot.save(plots, output)

def plot_violin(df, include_AA, group_col, x_label, cmap, background=False):
    """
    generates a list of either one or three violin plots depending on if AA mutation distributions are included

    df:             pd.DataFrame resulting from merged genotypes dataframes, then melted using the substitution count column(s)
    include_AA:     whether to include plots that describe AA mutation distribution in addition to NT
    """
    groups = list(df[group_col].unique())
    colors = get_colors(groups, cmap, background=background)
    df = df.rename(columns={'NT_substitutions_count':'NT','AA_substitutions_nonsynonymous_count':'AA'})

    # melt then reorder barcode groups back to original order
    value_cols = ['NT', 'AA'] if include_AA else ['NT']
    df = df.melt(id_vars=group_col, value_vars=value_cols, var_name='substitution_type', value_name='substitution_count')
    sorting_rule = {group: i for i,group in enumerate(groups)}
    df = df.sort_values(by=[group_col], key=lambda col: col.apply(lambda group: sorting_rule[group]))

    plots = []
    if include_AA:

        # hvplot messes with ordering for multi-categorical plots (sigh), so using hooks to explicitly define order
        def set_order(plot, element):
            labels = product([str(g) for g in groups], ['NT','AA'])
            plot.state.x_range.factors = [*labels]
        
        two_colors = [cmap_dict()[cmap][200], cmap_dict()[cmap][100]]
        both_plot = df.hvplot.violin(y='substitution_count', by=[group_col,'substitution_type'], c='substitution_type',
                        ylabel='substitutions per sequence', xlabel=x_label, legend=True, cmap=two_colors,
                        width=len(groups*100), height=500, title='substitutions distributions').opts(
                                hooks=[set_order], backend='bokeh')
        AA_plot = df[df['substitution_type']=='AA'].hvplot.violin(y='substitution_count', by=group_col, c=group_col,
                        ylabel='substitutions per sequence', xlabel=x_label, colormap=colors,
                        width=len(groups*50), height=500, legend=False, title='nonsynonymous AA substitutions distributions')
        plots.extend([both_plot,AA_plot])

    NT_plot = df[df['substitution_type']=='NT'].hvplot.violin(y='substitution_count', by=group_col, c=group_col,
                        ylabel='substitutions per sequence', xlabel=x_label, colormap=colors,
                        width=len(groups*50), height=500, legend=False, title='NT substitutions distributions')
    plots.append(NT_plot)

    return plots

if __name__ == '__main__':
    main(snakemake.input[0], snakemake.output[0], snakemake.params.group_col, snakemake.params.x_label, snakemake.params.export_SVG, snakemake.params.cmap, snakemake.params.background)