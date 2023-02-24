#
#  DESCRIPTION   : Script for maple pipeline. Uses nucleotide mutation data from mutation-stats.csv
#                   for designated timepoints to determine mutation rate and mutation spectrum.
#
#  AUTHOR(S)     : Gordon Rix
#

import os
import numpy as np
import pandas as pd
import bokeh
import re
import scipy
import math
from bokeh.layouts import column, gridplot
from bokeh.models import ColumnDataSource
from bokeh.palettes import Blues5
from bokeh.plotting import figure, output_file, save
from Bio import SeqIO
import holoviews as hv
import hvplot.pandas #noqa
from holoviews import opts
from scipy import stats
hv.extension('bokeh')

def nt_normal_dict(sequence):
    """ given a reference sequence, returns a dictionary of nucleotide:normalization factor key:value
        pairs such that multiplying this factor by (total mutations of type / total bases analyzed) will
        result in the mutation rate contributed by this type of mutation for a sequence with a 1:1:1:1 ratio
        of A:T:G:C
    """
    nt_normal_dict = {}
    for nt in 'ATGC':
        nt_normal_dict[nt] = len(sequence) / (4*sequence.count(nt))
    return nt_normal_dict

def trim_normalize_row(mutStatsRow, refSeq, mutTypes):
    """ given a row from mutation stats csv,
        returns the same row but trimmed to only include mutations per base, both total and of each type

        Args:
            mutStatsRow:        row from mut stats dataframe to be analyzed
            refSeq:             reference sequence used for mutation analysis
            mutTypes:           list of different types of mutations, used as column names
    """
    totalNTsAnalyzed = mutStatsRow['total_seqs'] * len(refSeq)
    mutList = []
    normalDict = nt_normal_dict(refSeq)
    for mut in mutTypes:
        wtNT = mut[0]
        mutList.append(normalDict[wtNT] * (mutStatsRow[mut]/totalNTsAnalyzed))

    # output total mutations per base as the sum of the normalized rates for all mutation types (may differ
    #   slightly from the rate as calculated by total mutations/total nts analyzed due to combination
    #   of mutation preference and A/T/G/C content of reference)
    outList = [sum(mutList)]
    outList.extend(mutList)

    return outList

def main():
    ### Asign variables from config file and inputs
    config = snakemake.config
    tag = snakemake.wildcards.tag
    mutStatsCSV = pd.read_csv(str(snakemake.input.mutStats))
    timepointsCSV = pd.read_csv(str(snakemake.input.timepoints), header=1, index_col=0).astype(str)
    topRow = [x for x in pd.read_csv(str(snakemake.input.timepoints)).columns if 'Unnamed: ' not in x]
    if len(topRow) > 0:
        timeUnit = topRow[0]
    else:
        timeUnit = 'generations'    # default label for time unit if none is provided in the first row of timepoints CSV
    backgroundBCgroup = config.get('background',False)
    refSeqfasta = config['runs'][tag]['reference']
    refSeq = str(list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq).upper()
    ###

    ### Make paths for outputs
    for outFile in snakemake.output:
        dir = '/'.join(outFile.split('/')[:-1])
        if not os.path.exists(dir):
            os.mkdir(dir)
    ###

    # generate list of all types of mutations to be used to filter dataframe to only include mutation type columns
    nts = 'ATGC'
    mutTypes = []
    for wt in nts:
        for mut in nts:
            if wt != mut:
                mutTypes.append(f'{wt}->{mut}')

    # dictionary of normalization factors for each nucleotide based on reference sequence
    nt_normal_dict = {}
    for nt in nts:
        nt_normal_dict[nt] = len(refSeq) / (4*refSeq.count(nt))

    # make a dictionary of tag:background rows from mutation stats csv file key:value pairs
    # Row is first trimmed for relevant stats and normalized by # of sequences and number of bases analyzed
    backgroundRows = {}
    if backgroundBCgroup:
        assert backgroundBCgroup in list(mutStatsCSV['barcode_group']), f'Provided barcode group for background subtraction, {backgroundBCgroup}, not present in {tag}_mutation-stats.csv. Demuxing for this barcode group may have failed.'
        for _, row in mutStatsCSV[mutStatsCSV['barcode_group']==backgroundBCgroup].iterrows():
            rowTag = row['tag']
            rowRefSeqfasta = config['runs'][rowTag]['reference']
            rowRefSeq = str(list(SeqIO.parse(rowRefSeqfasta, 'fasta'))[1].seq).upper()
            normTrimmedRow = trim_normalize_row(row, rowRefSeq, mutTypes)
            backgroundRows[rowTag] = normTrimmedRow

    allTimepointsDFrowList = [] # list to be populated with rows of mutations per sequence of each mutation type, one row for each timepoint
    allTimepointsDFcolumns = ['sample_label', 'replicate', timeUnit, 'per_base_all'] + ['per_base_'+mt for mt in mutTypes]
    allTimepointsDF = pd.DataFrame(allTimepointsDFrowList, columns=allTimepointsDFcolumns)

    allRatesDFrowList = [] # list to be populated with correlations between all types of mutations and timepoints, yielding a mutation rate as slope
    allRatesDFcolumns = ['sample_label', 'replicate', 'mut_type', 'wt_nt', 'mut_nt', 'rate', 'intercept', 'r_value', 'p_value', 'std_err']
    allRatesDF = pd.DataFrame(allRatesDFrowList, columns=allRatesDFcolumns)

    # loop through each row in the timepoints table, then loop through each column in the timepoints table,
    #   grab the total mutations as well as the total of each mutation type, divide by # of sequences,
    #   background subtract if background barcode group is given in config file, divide by # of bases per sequence
    #   (dependent on both sequence length and identity), yielding mutations per base analyzed of each type

    # loop through each unique sample label in timepoints, then loop through each timepoint for
    #   each of these samples (so that all replicates of a sample are handled at the same time)
    uniqueSamples = list(timepointsCSV.index.unique())
    for sampleLabel in uniqueSamples:

        replicateIndex = 0

        for _, row in timepointsCSV.loc[timepointsCSV.index==sampleLabel].iterrows():

            replicateIndex += 1
            sampleTimepointDFrowList = []

            for timepoint in timepointsCSV.columns:
                if row[timepoint]=='nan':
                    continue
                timepointTag, timepointBCgroup = row[timepoint].split('_')
                timepointRefSeqfasta = config['runs'][timepointTag]['reference']
                timepointRefSeq = str(list(SeqIO.parse(timepointRefSeqfasta, 'fasta'))[1].seq).upper()

                # grab row from mut stats corresponding to sample/barcode group
                timepointSeqsMutStatsRow = (mutStatsCSV.loc[(mutStatsCSV['tag']==timepointTag) & (mutStatsCSV['barcode_group']==timepointBCgroup)])

                if len(timepointSeqsMutStatsRow) == 1:
                    timepointSeqsMutStatsRow = (mutStatsCSV.loc[(mutStatsCSV['tag']==timepointTag) & (mutStatsCSV['barcode_group']==timepointBCgroup)])
                else:
                    print(f'Tag / barcodeGroup combination `{timepointTag}` / `{timepointBCgroup}` not present in mutStats CSV file. This timepoint for this sample will not be used to calculate mutation rates. Check demultiplexing definition and output for this sample.')
                    continue

                timepointSeqsMutStatsRow = timepointSeqsMutStatsRow.iloc[0]

                # calculate the substitutions per base analyzed for all types of substitutions individually and combined, normalized to # of sequences, and subtract sequencing background if given
                normTrimmedRow = trim_normalize_row(timepointSeqsMutStatsRow, timepointRefSeq, mutTypes)
                if backgroundBCgroup:
                    normTrimmedRow = list(np.array(normTrimmedRow) - np.array(backgroundRows[timepointTag]))
                sampleTimepointDFrowList.append([sampleLabel, replicateIndex, timepoint] + normTrimmedRow)
            
            if len(sampleTimepointDFrowList) < 2:
                print(f'Fewer than 2 timepoints identified for sample `{sampleLabel}` replicate `{replicateIndex}`. Mutation rate will not be calculated. Check demultiplexing definition and output for this sample.')
                continue

            sampleTimepointDF = pd.DataFrame(sampleTimepointDFrowList, columns=allTimepointsDFcolumns)
            sampleRatesDFrowList = []

            # calculate mutation rates for each type of mutation
            for mut_type in ['all'] + mutTypes:
                rate, intercept, r_value, p_value, std_err = scipy.stats.linregress(sampleTimepointDF[timeUnit].astype(float).values.tolist(), sampleTimepointDF['per_base_'+mut_type].astype(float).values.tolist())
                if mut_type == 'all':
                    wt, mut = 'all', 'all'
                else:
                    wt, mut = mut_type[0], mut_type[-1]
                sampleRatesDFrowList.append([sampleLabel, replicateIndex, mut_type, wt, mut, rate, intercept, r_value, p_value, std_err])

            allRatesDF = pd.concat([allRatesDF, pd.DataFrame(sampleRatesDFrowList, columns=allRatesDFcolumns)]).reset_index(drop=True)
            allTimepointsDF = pd.concat([allTimepointsDF, sampleTimepointDF]).reset_index(drop=True)

    # compute mean rates for replicate samples then remove negatives
    meanRatesDF = allRatesDF.groupby(['sample_label', 'mut_type', 'wt_nt', 'mut_nt'], sort=False)['rate'].describe().reset_index().rename(columns={'mean':'rate_mean', 'std':'rate_std'})
    meanRatesDF = meanRatesDF.drop(columns=meanRatesDF.columns[-5:])
    meanRatesDF['rate_mean'] = meanRatesDF['rate_mean'].clip(lower=10^-10) # convert negative values to 10^-10
    print(allRatesDF['rate'].sort_values()
    allRatesDF['rate'] = allRatesDF['rate'].clip(lower=10^-10)
    print(allRatesDF['rate'].sort_values()

    defaults = dict(height=400, tools=['hover'], fontsize={'title':16,'labels':14,'xticks':10,'yticks':10})
    boxwhisker_defaults = dict(box_fill_color='grey', box_line_width=1, whisker_line_width=1)
    hv.opts.defaults(hv.opts.BoxWhisker(**{**defaults, **boxwhisker_defaults}), hv.opts.HeatMap(**defaults))
    mutType_grouped_plot_list = []

    # plot that shows individual mutation rates, grouped together by mutation type (# of plots == # of mutation types including overall rate == 13)
    for mut_type in ['all'] + mutTypes:
        mtType_rate_DF = allRatesDF[allRatesDF['mut_type']==mut_type]
        boxPlot = hv.BoxWhisker(mtType_rate_DF, kdims='sample_label', vdims='rate').opts(
                    logy=True, xrotation=70, width=max(50*len(uniqueSamples), 150), # fails to render if too thin
                    xlabel='sample', outlier_alpha=0, # hide outliers because will show all points with Points
                    ylim=(0.000000001, 0.0005), ylabel=f'{mut_type} substitution rate')
        points = hv.Points(mtType_rate_DF, kdims='sample_label', vdims='rate').opts(
                    logy=True, color='black', alpha=0.7, jitter=0.2, size=6)
        mutType_grouped_plot_list.append(boxPlot*points)
        
    # plots that show individual rates, grouped together by sample (# of plots == # of samples)
    sample_grouped_plot_list = []
    heatmap_list = []
    for sample in uniqueSamples:

        # boxplot and points
        sample_rate_DF = allRatesDF[allRatesDF['sample_label']==sample]
        boxPlot = hv.BoxWhisker(sample_rate_DF, kdims=['wt_nt','mut_nt'], vdims='rate').opts(
                    logy=True, ylim=(0.000000001, 0.0005),
                    title=f'sample: {sample}', xlabel='mutation nucleotide\nwild type nucleotide',
                    width=650, ylabel=f'per base substitutions per {timeUnit}')
        points = hv.Points(mtType_rate_DF, kdims='sample_label', vdims='rate').opts(
                    logy=True, color='black', alpha=0.7, jitter=0.2, size=6)
        sample_grouped_plot_list.append(boxPlot*points)

        mean_rates_individual = meanRatesDF[(meanRatesDF['sample_label']==sample) & (meanRatesDF['wt_nt']!='all')
                                            ].sort_values(['wt_nt','mut_nt'], ascending=[True,False]) # sort in opposite order then flip yaxis to get same order for x and y axis
        heatmap = mean_rates_individual.hvplot.heatmap(x='wt_nt', y='mut_nt', C='rate_mean', by='sample_label',
                                                        flip_yaxis=True, width=480, title=f'sample: {sample}',
                                                        xlabel='wild type nucleotide', ylabel='mutation nucleotide'
                                                        ).opts(colorbar_opts={'title':f'substitutions per base per {timeUnit}'}, axiswise=True)
        heatmap_list.append(heatmap)

    # export rate plots and data
    hv.save(hv.Layout(mutType_grouped_plot_list).cols(1), snakemake.output.boxplot_mut_grouped, backend='bokeh')
    hv.save(hv.Layout(sample_grouped_plot_list).cols(1), snakemake.output.boxplot_plot_sample_grouped, backend='bokeh')
    hv.save(hv.Layout(heatmap_list).cols(1), snakemake.output.heatmap, backend='bokeh')
    allRatesDF.to_csv(snakemake.output.CSV_all_rates, index=False)

    # pivot mean rates, return to original order, then export
    meanRatesPivotDF = meanRatesDF.pivot(index='sample_label', columns='mut_type', values='rate_mean').reset_index()
    sortOrder = {value:position for position,value in enumerate(meanRatesDF['sample_label'].unique())}
    meanRatesPivotDF = meanRatesPivotDF.sort_values('sample_label', key=lambda x: x.apply(lambda y: sortOrder[y]))
    meanRatesPivotDF.to_csv(snakemake.output.CSV_summary, index=False)

if __name__ == '__main__':
    main()