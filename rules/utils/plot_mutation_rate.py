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
from common import export_svg_plots
from bokeh.layouts import column, gridplot
from bokeh.models import ColumnDataSource
from bokeh.palettes import Blues5
from bokeh.plotting import figure, output_file, save
from Bio import SeqIO
import holoviews as hv
import hvplot.pandas #noqa
from holoviews import opts
from scipy import stats
from common import cmap_dict
hv.extension('bokeh')

def nt_normal_dict(sequence):
    """ given a reference sequence, returns a dictionary of nucleotide:normalization factor key:value
        pairs such that multiplying this factor by (total mutations of type / total bases analyzed) will
        result in the mutations per base contributed by this type of mutation for a sequence with a 1:1:1:1 ratio
        of A:T:G:C
    """
    nt_normal_dict = {}
    for nt in 'ATGC':
        nt_normal_dict[nt] = len(sequence) / (4*sequence.count(nt))
    return nt_normal_dict

def trim_normalize_row(mutStatsRow, refSeq, mutTypes):
    """ given a row from mutation stats csv,
        returns the same row but trimmed to only include mutations and insertions per base, both total and of each type

        Args:
            mutStatsRow:        row from mut stats dataframe to be analyzed
            refSeq:             reference sequence used for mutation analysis
            mutTypes:           list of different types of mutations, used as column names
    """
    totalNTsAnalyzed = mutStatsRow['total_seqs'] * len(refSeq)

    # deletions can't be normalized to the sequence being looked at so calculate average insertion length per base first
    indel_list = []
    for indel in ['total_insertion_length', 'total_deletion_length']:
        if indel in mutTypes:
            mutTypes.remove(indel)
            indel_list.append(mutStatsRow[indel]/totalNTsAnalyzed)

    mutList = []
    normalDict = nt_normal_dict(refSeq)

    for mut in mutTypes:
        wtNT = mut[0]
        mutList.append(normalDict[wtNT] * (mutStatsRow[mut]/totalNTsAnalyzed))

    # output total mutations per base as the sum of the normalized rates for all mutation types (may differ
    #   slightly from the rate as calculated by total mutations/total nts analyzed due to combination
    #   of mutation preference and A/T/G/C content of reference)
    outList = [sum(mutList)]
    outList.extend(mutList+indel_list)

    return outList

def horizontal_lines(start, end):
    """
    Generates a holoviews overlay of horizontal lines at powers of 10 in the provided range
        (10^n where n = range(start,end)

    Args:
        start:      start of range of exponents for horizontal lines
        end:        end of range of exponents for horizontal lines

    Returns:
        holoviews overlay plot
    
    """
    # Ensure we have a range function for numpy
    exponent_range = [10**n for n in range(start,end)]
    lines = []
    for y in exponent_range:
        lines.append( hv.HLine(y=y).opts(color="lightgrey", line_dash="dashed") )
    return hv.Overlay(lines)

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
        timeUnit = 'unit time'    # default label for time unit if none is provided in the first row of timepoints CSV
    backgroundBCgroup = config.get('background',False)
    refSeqfasta = config['runs'][tag]['reference']
    refSeq = str(list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq).upper()
    cmap = cmap_dict()[snakemake.params.cmap]
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
            normTrimmedRow = trim_normalize_row(row, rowRefSeq, mutTypes+['total_insertion_length', 'total_deletion_length'])
            backgroundRows[rowTag] = normTrimmedRow

    allTimepointsDFcolumns = ['sample_label', 'replicate', timeUnit, 'per_base_all'] + ['per_base_'+mt for mt in mutTypes] + ['per_base_insertion', 'per_base_deletion']
    allTimepointsDF = pd.DataFrame([], columns=allTimepointsDFcolumns)

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
                normTrimmedRow = trim_normalize_row(timepointSeqsMutStatsRow, timepointRefSeq, mutTypes+['total_insertion_length', 'total_deletion_length'])
                if backgroundBCgroup:
                    normTrimmedRow = list(np.array(normTrimmedRow) - np.array(backgroundRows[timepointTag]))
                sampleTimepointDFrowList.append([sampleLabel, replicateIndex, timepoint] + normTrimmedRow)
            
            if len(sampleTimepointDFrowList) < 2:
                print(f'Fewer than 2 timepoints identified for sample `{sampleLabel}` replicate `{replicateIndex}`. Mutation rate will not be calculated. Check demultiplexing definition and output for this sample.')
                continue

            sampleTimepointDF = pd.DataFrame(sampleTimepointDFrowList, columns=allTimepointsDFcolumns)
            sampleRatesDFrowList = []

            # calculate mutation rates for each type of mutation
            for mut_type in ['all'] + mutTypes + ['insertion', 'deletion']:
                rate, intercept, r_value, p_value, std_err = scipy.stats.linregress(sampleTimepointDF[timeUnit].astype(float).values.tolist(), sampleTimepointDF['per_base_'+mut_type].astype(float).values.tolist())
                if mut_type not in mutTypes:
                    wt, mut = mut_type[:3], mut_type[:3] # all, ins, del
                else:
                    wt, mut = mut_type[0], mut_type[-1]
                sampleRatesDFrowList.append([sampleLabel, replicateIndex, mut_type, wt, mut, rate, intercept, r_value, p_value, std_err] + sampleTimepointDF['per_base_'+mut_type].astype(float).values.tolist())

            # make new DF that includes all normalized mutations per base
            allRatesDF = pd.concat([allRatesDF, pd.DataFrame(sampleRatesDFrowList, columns=allRatesDFcolumns + sampleTimepointDF[timeUnit].astype(float).values.tolist())]).reset_index(drop=True)
            allTimepointsDF = pd.concat([allTimepointsDF, sampleTimepointDF]).reset_index(drop=True)

    # remove negatives then compute mean rates for replicate samples 
    allRatesDF['rate'] = allRatesDF['rate'].clip(lower=10**-10)
    meanRatesDF = allRatesDF.groupby(['sample_label', 'mut_type', 'wt_nt', 'mut_nt'], sort=False)['rate'].describe().reset_index().rename(columns={'mean':'rate_mean', 'std':'rate_std'})
    meanRatesDF = meanRatesDF.drop(columns=meanRatesDF.columns[-5:])
    meanRatesDF['rate_mean'] = meanRatesDF['rate_mean'].clip(lower=10**-8) # convert means <= 10^-8 to nan bc they are not reliable
    meanRatesDF.loc[meanRatesDF['rate_mean'] == 10**-8, 'rate_mean'] = np.nan

    def hook(plot, element):
        plot.output_backend = 'svg'

    defaults = dict(height=400, tools=['hover'], fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}, hooks=[hook])
    boxwhisker_defaults = dict(box_fill_color='grey', box_line_width=1, whisker_line_width=1)
    hv.opts.defaults(hv.opts.BoxWhisker(**{**defaults, **boxwhisker_defaults}), hv.opts.HeatMap(**defaults))
    mutType_grouped_plot_list = []

    # plot that shows individual mutation rates, grouped together by mutation type (# of plots == # of mutation types including overall rate and in/dels == 15)
    for mut_type in ['all'] + mutTypes + ['insertion', 'deletion']:
        mtType_rate_DF = allRatesDF[allRatesDF['mut_type']==mut_type]
        if mut_type in ['all']+mutTypes:
            ylabel = f'per base {mut_type} substitution rate'
        else:
            ylabel = f'per base {mut_type} rate'
        boxPlot = hv.BoxWhisker(mtType_rate_DF, kdims='sample_label', vdims='rate').opts(
                    logy=True, xrotation=70, width=max(50*len(uniqueSamples), 150), # fails to render if too thin
                    xlabel='sample', outlier_alpha=0, # hide outliers because will show all points with Points
                    ylim=(0.00000001, 0.0005), ylabel=ylabel)
        points = hv.Points(mtType_rate_DF[['sample_label', 'rate']]).opts(
                    logy=True, color='black', alpha=0.7, jitter=0.2, size=6, tools=['hover'])
        # 10^n markers
        h_lines = horizontal_lines(-11, -1)

        mutType_grouped_plot_list.append(h_lines*boxPlot*points)
        
    # plots that show individual rates, grouped together by sample (# of plots == # of samples)
    sample_grouped_plot_list = []
    heatmap_list = []
    for sample in uniqueSamples:

        # no points plot because points can't do multi category x axis
        sample_rate_DF = allRatesDF[allRatesDF['sample_label']==sample]
        boxPlot = hv.BoxWhisker(sample_rate_DF, kdims=['wt_nt','mut_nt'], vdims='rate').opts(
                    logy=True, ylim=(0.00000001, 0.0005),
                    title=f'sample: {sample}', xlabel='mutation nucleotide\nwild type nucleotide',
                    width=650, ylabel=f'per base substitutions per {timeUnit}', tools=['hover'])
        # points = hv.Points(sample_rate_DF, kdims=['wt_nt','mut_nt'], vdims='rate').opts(
        #             logy=True, color='black', alpha=0.7, jitter=0.2, size=6)
        # 10^n markers (n=integer)
        h_lines = horizontal_lines(-11, -1)
        sample_grouped_plot_list.append(h_lines*boxPlot)

        # plot individual substitution rates as heatmap
        mean_rates_individual = meanRatesDF[(meanRatesDF['sample_label']==sample) & (~meanRatesDF['wt_nt'].isin(['all','ins','del']))
                                            ].sort_values(['wt_nt','mut_nt'], ascending=[True,False]) # sort in opposite order then flip yaxis to get same order for x and y axis
        mean_rates_individual['rate_mean'] = np.log10(mean_rates_individual['rate_mean']) # log10 transform, but note that only the integer tick marks on the color bar will be valid now
        heatmap = mean_rates_individual.hvplot.heatmap(x='wt_nt', y='mut_nt', C='rate_mean', by='sample_label',
                                                        flip_yaxis=True, width=480, title=f'sample: {sample}', cmap=cmap,
                                                        xlabel='wild type nucleotide', ylabel='mutation nucleotide', clim=(-8, np.log10(0.00015),
                                                        ).opts(colorbar_opts={'title':f'log10(s.p.b per {timeUnit})'}, axiswise=False)
        heatmap_list.append(heatmap)

    # export rate plots and data
    hv.save(hv.Layout(mutType_grouped_plot_list).cols(1), snakemake.output.boxplot_mut_grouped, backend='bokeh')
    hv.save(hv.Layout(sample_grouped_plot_list).cols(1), snakemake.output.boxplot_plot_sample_grouped, backend='bokeh')
    hv.save(hv.Layout(heatmap_list).cols(1), snakemake.output.heatmap, backend='bokeh')

    if snakemake.params.export_SVG:
        export_svg_plots(mutType_grouped_plot_list, snakemake.output.boxplot_mut_grouped, labels=['all']+mutTypes, export=snakemake.params.export_SVG)
        export_svg_plots(sample_grouped_plot_list, snakemake.output.boxplot_plot_sample_grouped, labels=uniqueSamples, export=snakemake.params.export_SVG)
        export_svg_plots(heatmap_list, snakemake.output.heatmap, labels=uniqueSamples, export=snakemake.params.export_SVG) # note: colorbars can't be exported to SVG

    allRatesDF.to_csv(snakemake.output.CSV_all_rates, index=False)

    # pivot mean rates, return to original order, then export
    meanRatesPivotDF = meanRatesDF.pivot(index='sample_label', columns='mut_type', values='rate_mean').reset_index()
    sortOrder = {value:position for position,value in enumerate(meanRatesDF['sample_label'].unique())}
    meanRatesPivotDF = meanRatesPivotDF.sort_values('sample_label', key=lambda x: x.apply(lambda y: sortOrder[y]))
    meanRatesPivotDF.to_csv(snakemake.output.CSV_summary, index=False)

if __name__ == '__main__':
    main()