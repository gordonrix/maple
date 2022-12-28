#
#  DESCRIPTION   : Script for maple pipeline. Generates bar plots for hamming distance
#                       distributions
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import holoviews as hv
import hvplot.pandas
import colorcet as cc
hv.extension('bokeh')

def main(input, output, labels, title, legendLabel, background, raw):

    colormap = cc.blues
    colorConstant = (len(colormap)-1) / len(input)
    colors = []

    for i, CSV in enumerate(input):
        _, barcodes, _ = CSV.split('/')[-1].split('_')
        color = colormap[int(i*colorConstant)]
        if background:
            if barcodes == background:
                color = 'grey'
        colors.append(color)

    inputDFs = [pd.read_csv(CSV, index_col=False) for CSV in input]
    plots = [ plot_cumsum(inputDFs, labels, colors, title, legendLabel) ] + [
              plot_dist(df, title=f"{legendLabel}: {label}", raw=raw) for df, label in zip(inputDFs, labels) ]
    plots = hv.Layout(plots).cols(1)
    hvplot.save(plots, output)

def plot_cumsum(DFlist, labels, colors, title, legendLabel):
    """
    generates a distribution hv.Curve plot showing cumulative sums for all input files as hv.Curve on an hv.Overlay plot
    
    DFlist:         list of pd.DataFrames that contain cumulative sums in the format generated by the common.dist_to_DF function
    labels:         list of strings, names to use to label each group for cumulative distribution hv.Overlay
    colors:         list of strings, color values index matched to DFlist to be used to color each line
    title:          string, label for all samples to include in the title
    legendLabel:    string, label to add to the legend (only appears in hover)
    """

    plotList = []
    for dfOG, label, color in zip(DFlist, labels, colors):

        df = dfOG.copy()    # prevent adding a column to the global dataframe
        x, total, proportion, cumsum = df.columns
        y = total.split('total ')[1]
        df[legendLabel] = label

        plotList.append( hv.Curve(df, kdims=[x],
                            vdims = [cumsum, proportion, total] + [legendLabel], label=label).opts(
                                        hv.opts.Curve(tools=['hover'], color=color)) )

    plot = hv.Overlay(plotList).opts(show_legend=True, legend_position='right', width=800, height=600,
                            title=f"{x} among {y}, cumulative distribution\nsample: {title}",
                            # legend_title = legendLabel, # This doesn't work for some reason, but a legend title would be nice
                            xlabel=x, ylabel='cumulative frequency',
                            fontsize={'title':16,'labels':14,'xticks':12,'yticks':12})

    return plot


def plot_dist(df, color='grey', title='', raw=False):
    """
    generates a histogram with bin size of 0 for a distribution input file

    df:     pd.DataFrame, distribution data in the format generated by the common.dist_to_DF function
    color:  string, color to use for bars, all will be the same color
    title:  string, title that describes the individual sample
    raw:    bool, determines if y axis will show raw counts or proportion of the total
                default behavior is to show proportion of the total, not raw counts
    """

    x, total, proportion, cumsum = df.columns
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
    plot = hv.Histogram(hvDF).opts(hv.opts.Histogram(tools=['hover'], color=color, width=800, height=600,
                                    title=f"{x} among {y}\n{title}",
                                    ylabel=ylabel,
                                    fontsize={'title':16,'labels':14,'xticks':12,'yticks':12}))

    return plot

if __name__ == '__main__':
    main(snakemake.input, snakemake.output.plot, snakemake.params.labels, snakemake.params.title, snakemake.params.legend_label, snakemake.params.background, snakemake.params.raw)