#
#  DESCRIPTION   : Script for maple pipeline. Uses nucleotide mutation data from mutation-stats.csv
#                   to determine mutation spectrum.
#
#  AUTHOR(S)     : Gordon Rix
#

import os
import numpy as np
import pandas as pd
import re
from bokeh.layouts import column
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource, FactorRange,
                          HoverTool, Legend, LinearColorMapper,
                          PrintfTickFormatter)
from bokeh.palettes import Blues5
from bokeh.plotting import figure, output_file, save, show
from Bio import SeqIO

def normalized_spectrum_df(mutStatsRow, backgroundRow, backgroundBool, mutTypes, nt_normal_dict):
    """ returns a dataframe of length 1 with a column for each type of substitution designating the proportion
    of mutations of that type of mutation. This value is normalized based on the frequency of appearance in the
    reference sequence, and is background subtracted if a background is provided
    Args:
        mutStatsRow:        row from mut stats dataframe to be analyzed
        backgroundRow:      row from mut stats dataframe for background sequences
        backgroundBool:     boolean used to determine if background subtraction should be done
        mutTypes:           list of different types of mutations
        nt_normal_dict:     dictionary of nucleotide:normalization factor key:value pairs
    """

    spectrumNormalizedDict = {}

    for mutType in mutTypes:
        wtNT = mutType[0]
        if mutStatsRow['total_NT_mutations'] == 0:  # set all mutation types to 0 if no mutations are observed to avoid divide by 0
            spectrumNormalizedDict[mutType] = [ 0 ]
        elif backgroundBool:
            spectrumNormalizedDict[mutType] = [ max(0, nt_normal_dict[wtNT] * ( (mutStatsRow[mutType] / mutStatsRow['total_seqs']) - (backgroundRow[mutType] / backgroundRow['total_seqs']) )) ]    # reference sequence normalization factor * ((number of type of mutation for row / total number of sequences for row) - (number of type of mutation for background / total number of sequences for background)) if negative, will be 0
        else:
            spectrumNormalizedDict[mutType] = [ nt_normal_dict[wtNT] * (mutStatsRow[mutType] / mutStatsRow['total_NT_mutations']) ]    # reference sequence normalization factor * (number of type of mutation for row / total number of mutations for row)

    outDF = pd.DataFrame(spectrumNormalizedDict).reset_index(drop=True)
    finalNorm = 1 / outDF.to_numpy().sum()
    outDF = outDF * finalNorm

    return outDF

def main():
    ### Asign variables from config file and inputs
    config = snakemake.config
    tag = snakemake.wildcards.tag
    mutStats = pd.read_csv(str(snakemake.input))
    background = config.get('background',False)
    refSeqfasta = config['runs'][tag]['reference']
    refSeq = str(list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq).upper()
    ###

    ### Output variables
    plotOut = snakemake.output[0]
    csvOutDir = 'mutation_data'
    if not os.path.exists(csvOutDir):
        os.mkdir(csvOutDir)
    ###

    # generate list of all types of mutations to be used to filter dataframe to only include mutation type columns
    nts = 'ATGC'
    mutTypes = []
    for wt in nts:
        for mut in nts:
            if wt != mut:
                if config.get('uniques_only', False):   col = f'{wt}->{mut}_unique'
                else:                                   col = f'{wt}->{mut}'
                mutTypes.append(col)

    # dictionary of normalization factors for each nucleotide based on reference sequence
    nt_normal_dict = {}
    for nt in nts:
        nt_normal_dict[nt] = len(refSeq) / (4*refSeq.count(nt))

    if background:
        assert background in list(mutStats['barcode_group']), f'Provided barcode group for background subtraction, {background}, not present in {tag}_mutation-stats.csv'
        for _, row in mutStats[mutStats['barcode_group']==background].iterrows():
            backgroundRow = row
        backgroundBool = True
    else:
        backgroundRow = None
        backgroundBool = False

    output_file(plotOut)
    plotList = []
        
    first = True
    # generate spectrum DF for each barcode in the sample, reformat for plotting, then plot as column
    for _, row in mutStats.iterrows():

        if row['total_seqs'] == 0: continue

        if backgroundBool:
            if backgroundRow['barcode_group'] == row['barcode_group']:  # don't background subtract for background sample
                backgroundBool = False

        bcGroupSpectrum = normalized_spectrum_df(row, backgroundRow, backgroundBool, mutTypes, nt_normal_dict)

        bcGroupSpectrumMelted = pd.melt(bcGroupSpectrum, value_vars=bcGroupSpectrum)
        bcGroupSpectrumMelted['wtNT'] = bcGroupSpectrumMelted.apply(lambda row:
            row['variable'][0], axis=1)
        bcGroupSpectrumMelted['mutNT'] = bcGroupSpectrumMelted.apply(lambda row:
            row['variable'][3], axis=1)
        
        bcGroupSpectrumPivot = bcGroupSpectrumMelted.pivot(index='wtNT', columns='mutNT', values = 'value')

        dataOutName = os.path.join(csvOutDir, tag, row.barcode_group, f'{tag}_{row.barcode_group}_mutSpectrum.csv')
        renameDict = {}
        for nt in list('ATGC'): renameDict[nt] = 'mut: '+nt
        bcGroupSpectrumPivot.rename(columns=renameDict).to_csv(dataOutName)

        bcGroupSpectrumPivot = bcGroupSpectrumPivot.fillna(0.0) # replace NaNs with 0
        bcGroupSpectrumPivot = bcGroupSpectrumPivot.clip(lower=0) # convert negative values to 0

        TOOLTIPS = [('mutation type', '@wtNT'+'->'+'$name'), ('proportion of mutations', '@$name')]
        plotTitle = f"{tag}_{row['barcode_group']}"
        spectrumPlot = figure(title=plotTitle, plot_height=400, plot_width=400, x_range=list('ATGC'), tooltips=TOOLTIPS)

        if first: # legend only for first plot
            spectrumPlot.vbar_stack(list('ATGC'), x='wtNT', width = 0.9, source=ColumnDataSource(bcGroupSpectrumPivot),
                color = Blues5[:4], legend_label=list('ATGC'))
            spectrumPlot.legend.title = 'mutation'
            first = False
        else:
            spectrumPlot.vbar_stack(list('ATGC'), x='wtNT', width = 0.9, source=ColumnDataSource(bcGroupSpectrumPivot),
                color = Blues5[:4])
        
        spectrumPlot.yaxis.axis_label = 'normalized substitution frequency'

        plotList.append(spectrumPlot)

    if len(plotList)>0:
        plotList[-1].xaxis.axis_label = 'wild type nucleotide' # label x axis for bottom plot

    save(column(plotList))

if __name__ == '__main__':
    main()