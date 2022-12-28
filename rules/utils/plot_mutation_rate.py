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
from bokeh.layouts import column, gridplot
from bokeh.models import ColumnDataSource
from bokeh.palettes import Blues5
from bokeh.plotting import figure, output_file, save
from Bio import SeqIO
import holoviews as hv
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

    ### Output variables
    ratePlotOut = snakemake.output.rate
    rateCSVout = snakemake.output.rateCSV
    spectrumPlotOut = snakemake.output.spectrum
    spectrumCSVout = snakemake.output.spectrumCSV
    for outFile in [ratePlotOut, rateCSVout, spectrumPlotOut, spectrumCSVout]:
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

    allTimepointsDFrowList = [] # list to be populated with one row as a list for each timepoint
    allTimepointsDFcolumns = ['sample_label', 'replicate', timeUnit, 'per_base_all'] + ['per_base_'+mt for mt in mutTypes]
    allTimepointsDF = pd.DataFrame(allTimepointsDFrowList, columns=allTimepointsDFcolumns)

    allRatesDFrowList = [] # list to be populated with correlations between all types of mutations and timepoints, yielding a mutation rate as slope
    allRatesDFcolumns = ['sample_label', 'replicate', 'rate_type', 'rate', 'intercept', 'r_value', 'p_value', 'std_err']
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
            for mt in ['per_base_all'] + ['per_base_'+mt for mt in mutTypes]:
                rate, intercept, r_value, p_value, std_err = scipy.stats.linregress(sampleTimepointDF[timeUnit].astype(float).values.tolist(), sampleTimepointDF[mt].astype(float).values.tolist())
                sampleRatesDFrowList.append([sampleLabel, replicateIndex, mt+f'_substitutions_per_{timeUnit}', rate, intercept, r_value, p_value, std_err])

            allRatesDF = pd.concat([allRatesDF, pd.DataFrame(sampleRatesDFrowList, columns=allRatesDFcolumns)]).reset_index(drop=True)
            allTimepointsDF = pd.concat([allTimepointsDF, sampleTimepointDF]).reset_index(drop=True)

    # compute mean rates for replicate samples then remove negatives
    meanRatesDF = allRatesDF.groupby(['sample_label', 'rate_type'])['rate'].mean().reset_index()
    meanRatesDF['rate'] = meanRatesDF['rate'].clip(lower=10^-9) # convert negative values to 10^-9

    defaults = dict(width=100*len(uniqueSamples), xrotation=70, height=400, tools=['tap', 'hover', 'box_select'])
    hv.opts.defaults(hv.opts.Bars(**defaults), hv.opts.Points(**defaults))
    ratePlotList = []

    for mt in ['all'] + mutTypes:
        # plot mean as bar
        mtTypeMeanRateDF = meanRatesDF[meanRatesDF['rate_type']==f'per_base_{mt}_substitutions_per_{timeUnit}']
        meanPlot = hv.Bars(mtTypeMeanRateDF[['sample_label','rate']])
        meanPlot.opts(logy=True, color='grey')

        # plot individual replicates as points
        mtTypeRepsRateDF = allRatesDF[allRatesDF['rate_type']==f'per_base_{mt}_substitutions_per_{timeUnit}']
        repsPlot = hv.Points(mtTypeRepsRateDF[['sample_label','rate']])
        repsPlot.opts(logy=True, color='black', alpha=0.7, jitter=0.2, size=6)

        plot = meanPlot * repsPlot

        # log plots in holoviews currently seem to be quite buggy, so I am using hardcoded y axis bounds.
        #   Might be able to make this better adapted to data if holoviews/bokeh fix this. submitted a bug report.
        plot.opts(ylim=(0.0000009, 0.0005))
        plot.opts(ylabel=f'per_base_{mt}_substitutions_per_{timeUnit}')
        ratePlotList.append(plot)
        
    # export rate plot and data for replicates
    if len(ratePlotList) > 0:
        hv.save(hv.Layout(ratePlotList).cols(1), ratePlotOut, backend='bokeh')
    allRatesDF.to_csv(rateCSVout)
    
    # DF that will convert the different absolute mutation rates into relative rates
    relativeSpectrumDF = allRatesDF
    # make new column that includes both sample name and replicate
    relativeSpectrumDF['sample_replicate'] = relativeSpectrumDF.apply(lambda row:
        str(row['sample_label'])+'_'+str(row['replicate']), axis=1)
    
    # rename column and column values to serve as better column titles after pivot, e.g. A->T
    relativeSpectrumDF['mutation_type'] = relativeSpectrumDF.apply(lambda row:
        row['rate_type'].split('_')[2], axis=1)
    
    relativeSpectrumDF = relativeSpectrumDF.pivot(index='sample_replicate', columns='mutation_type', values='rate')
    # zero out negative rates, compute new total rates, compute relative rates, remove 'all' column
    relativeSpectrumDF[relativeSpectrumDF<0] = 0
    relativeSpectrumDF['all'] = relativeSpectrumDF.apply(lambda row:
        sum(row[mutTypes]), axis=1)
    for mt in mutTypes:
        relativeSpectrumDF[mt] = relativeSpectrumDF.apply(lambda row:
            row[mt]/row['all'], axis=1)

    relativeSpectrumDF.drop(columns='all', inplace=True)

    # melt to manipulate column names
    relativeSpectrumDF = relativeSpectrumDF.melt(value_vars=relativeSpectrumDF.columns, ignore_index=False)
    relativeSpectrumDF['wtNT'] = relativeSpectrumDF.apply(lambda row:
        row['mutation_type'][0], axis=1)
    relativeSpectrumDF['mutNT'] = relativeSpectrumDF.apply(lambda row:
        row['mutation_type'][3], axis=1)

    # dictionary for renaming columns
    renameDict = {}
    for nt in list('ATGC'): renameDict[nt] = 'mut: '+nt
    spectrumPlotList = []

    for sample_replicate in relativeSpectrumDF.index.unique():
        plotDF = relativeSpectrumDF[relativeSpectrumDF.index==sample_replicate]
        plotDF = plotDF.pivot(index='wtNT', columns='mutNT', values = 'value')
        plotDF.rename(columns=renameDict)
        plotDF = plotDF.fillna(0.0) # replace NaNs with 0

        TOOLTIPS = [('mutation type', '@wtNT'+'->'+'$name'), ('proportion of mutations', '@$name')]
        spectrumPlot = figure(title=sample_replicate, plot_height=400, plot_width=460, x_range=list('ATGC'), tooltips=TOOLTIPS)

        spectrumPlot.vbar_stack(list('ATGC'), x='wtNT', width = 0.9, source=ColumnDataSource(plotDF),
            color = Blues5[:4], legend_label=list('ATGC'))
        spectrumPlot.legend.title = 'mutation'
        spectrumPlot.add_layout(spectrumPlot.legend[0], 'right')

        spectrumPlot.yaxis.axis_label = 'normalized substitution frequency'

        spectrumPlotList.append(spectrumPlot)

        spectrumPlotList[-1].xaxis.axis_label = 'wild type nucleotide' # label x axis for bottom plot

    output_file(spectrumPlotOut)
    if len(spectrumPlotList) > 0:
        save(column(spectrumPlotList))
    relativeSpectrumDF.to_csv(spectrumCSVout)

if __name__ == '__main__':
    main()