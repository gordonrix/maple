#
#  DESCRIPTION   : Script for maple pipeline. Generates bar plots that describe the demultiplexing outcomes,
#                   including some extensive plottting for noSplit barcodes
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import colorcet as cc
import holoviews as hv
from natsort import index_natsorted
from common import dist_to_DF
from plot_distribution import plot_cumsum, plot_dist
import hvplot.pandas

def main(input, output, barcodeInfoDict, barcodeGroupsDict):
    demuxDF = pd.read_csv(input)
    assert len(demuxDF['tag'].unique())==1   # this script uses demux stats output for only a single tag
    demuxDF = demuxDF.drop(columns='tag')

    color = 'grey'
    colormap = cc.blues

    noSplitBCtypes = [bcType for bcType in barcodeInfoDict.keys() if barcodeInfoDict[bcType].get('noSplit', False)]
    groupBCtypes = [bcType for bcType in barcodeInfoDict.keys() if bcType not in noSplitBCtypes]

    bcGroups = list(barcodeGroupsDict.keys()) if barcodeGroupsDict else []
    
    # group according to output file names
    aggs = {'output_file_barcodes':'first', 'demuxed_count':'first'}
    aggs.update({bcType:'first' for bcType in groupBCtypes})
    dfGrouped = demuxDF.groupby('output_file_barcodes', as_index=False).agg(aggs).rename(columns={'demuxed_count':'total_barcodes_count'})

    # split into groups with user-provided names vs not, sort by user-provided names and natural sort, respectively, then combine
    dfGroupedNames = dfGrouped[ dfGrouped['output_file_barcodes'].isin(bcGroups) ]
    dfGroupedNames = dfGroupedNames.sort_values( by='output_file_barcodes', key= lambda col: [bcGroups.index(x) for x in col] )
    dfGroupedNoNames = dfGrouped[ ~dfGrouped['output_file_barcodes'].isin(bcGroups)]
    dfGroupedNoNames = dfGroupedNoNames.sort_values( by='output_file_barcodes', key= lambda x: np.argsort(index_natsorted(dfGroupedNoNames['output_file_barcodes'])) )

    dfGrouped = pd.concat([dfGroupedNames, dfGroupedNoNames])

    plot = dfGrouped.hvplot.bar(x='output_file_barcodes', y='total_barcodes_count', hover_cols=['fwd','rvs'],
                                color=color, title="total demultiplex counts for each group", rot=70,
                                fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}, width=800, height=600)

    if noSplitBCtypes:
        
        # find all rows in which any noSplit barcodes failed identification, then group to get the total number, then add this as a column to dfGrouped
        query = " | ".join([x + " == 'fail'" for x in noSplitBCtypes])
        dfFailedGrouped = demuxDF.query(query).groupby('output_file_barcodes', as_index=False).agg({'barcodes_count':'sum'}).rename(columns={'barcodes_count':'failed_noSplit_barcodes_count'})
        dfGrouped = dfGrouped.merge(dfFailedGrouped, how='left', on='output_file_barcodes')
        dfGrouped['failed_noSplit_barcodes_count'] = dfGrouped['failed_noSplit_barcodes_count'].fillna(0).astype(int)
        dfGrouped['success_noSplit_barcodes_count'] = dfGrouped['total_barcodes_count'] - dfGrouped['failed_noSplit_barcodes_count']

        # plot number of sequences in each group with one or more failed noSplit barcode identification
        plot = plot + dfGrouped.hvplot.bar(x='output_file_barcodes', y='failed_noSplit_barcodes_count', hover_cols=['fwd','rvs'],
                                            color=color, title="failed noSplit demultiplex counts for each group", rot=70,
                                            fontsize={'title':16,'labels':14,'xticks':12,'yticks':12}, width=800, height=600)
        
        # plot number of sequences in each group with successful identification of all noSplit barcodes
        plot = plot + dfGrouped.hvplot.bar(x='output_file_barcodes', y='success_noSplit_barcodes_count', hover_cols=['fwd','rvs'],
                                            color=color, title="successful noSplit demultiplex counts for each group", rot=70,
                                            fontsize={'title':16,'labels':14,'xticks':12,'yticks':12}, width=800, height=600)

        # get distributions of noSplit barcode counts for each barcode group, then plot
        distList = []
        labels = []
        colors = []
        colorConstant = (len(colormap)-1) / len(dfGrouped['output_file_barcodes'].unique())

        for i, output_file_barcode in enumerate(dfGrouped['output_file_barcodes']):
            subset = demuxDF[demuxDF['output_file_barcodes']==output_file_barcode]
            bincounts = np.bincount(subset['barcodes_count'])
            distList.append( dist_to_DF(np.trim_zeros(bincounts,'b'), 'noSplit demux count', 'unique barcode groups') )
            labels.append(output_file_barcode)
            colors.append( colormap[int(i*colorConstant)] )

        legendLabel = 'output file group'
        distPlots = [ plot_cumsum( distList, labels, colors, 'noSplit demux', legendLabel ) ] + [
                 plot_dist(df, title=f"{legendLabel}: {label}", raw=True) for df, label in zip(distList, labels) ]
        distPlots = hv.Layout(distPlots)

        plot = plot + distPlots
        plot = plot.cols(1)

    hvplot.save(plot, output)

if __name__ == '__main__':
    main(snakemake.input.CSV, snakemake.output.plot, snakemake.params.barcodeInfo, snakemake.params.barcodeGroups)