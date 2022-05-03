""" script for nanoMACE pipeline

Uses nucleotide mutation data from mutation-stats.csv file to determine normalized (and background subtracted if specified)
mutation spectrum. Outputs mutation spectrum as bokeh .html plot, and the normalized data as a .csv (silently because
these outputs have different wildcards, which snakemake doesn't tolerate
"""


import os
import numpy as np
import pandas as pd
import seaborn as sns
import re
from bokeh.layouts import column
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource, FactorRange,
                          HoverTool, Legend, LinearColorMapper,
                          PrintfTickFormatter)
from bokeh.palettes import Blues5
from bokeh.plotting import figure, output_file, save, show
from Bio import SeqIO
import holoviews as hv
from holoviews import opts
from holoviews.streams import Selection1D
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

    # calculate mutations per base as the sum of the normalized rates for all mutation types (may differ from the rate
    # as calculated by total mutations/total nts analyzed due to mutation preferences and A/T/G/C content of reference
    outList = [sum(mutList)]
    outList.extend(mutList)

    return outList

def main():
    ### Asign variables from config file and inputs
    config = snakemake.config
    tag = snakemake.wildcards.tag
    mutStatsCSV = pd.read_csv(str(snakemake.input.mutStats))
    timepointsCSV = pd.read_csv(str(snakemake.input.timepoints), header=1, index_col=0)
    topRow = [x for x in pd.read_csv(str(snakemake.input.timepoints)).columns if 'Unnamed: ' not in x]
    if len(topRow) > 0:
        timeUnit = topRow[0]
    else:
        timeUnit = 'generations'    # use as label for time unit if none is provided in the first row of timepoints CSV
    backgroundBCgroup, backgroundBool = (config['background'], True) if 'background' in config else (None, False)
    refSeqfasta = config['runs'][tag]['reference']
    refSeq = str(list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq).upper()
    ###

    ### Output variables
    ratePlotOut = snakemake.output.rate
    # spectrumPlotOut = snakemake.output.spectrum
    # rateCSVout = snakemake.output.rateCSV
    # for outFile in [ratePlotOut, spectrumPlotOut, rateCSVout]:
    #     dir = '/'.join(outFile.split('/')[:-1])
    #     if not os.path.exists(dir):
    #         os.mkdir(dir)
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
    if backgroundBool:
        assert backgroundBCgroup in list(mutStatsCSV['barcode_group']), f'Provided barcode group for background subtraction, {backgroundBCgroup}, not present in {tag}_mutation-stats.csv. Demuxing for this barcode group may have failed.'
        for _, row in mutStatsCSV[mutStatsCSV['barcode_group']==backgroundBCgroup].iterrows():
            rowTag = row['tag']
            rowRefSeqfasta = config['runs'][rowTag]['reference']
            rowRefSeq = str(list(SeqIO.parse(rowRefSeqfasta, 'fasta'))[1].seq).upper()
            normTrimmedRow = trim_normalize_row(row, rowRefSeq, mutTypes)
            backgroundRows[rowTag] = normTrimmedRow

    sampleTimepointDFrowList = [] # list to be populated with one row as a list for each timepoint
    sampleTimepointDFcolumns = ['sampleLabel', 'replicate', 'generations', 'mutations_per_base'] + ['per_base_'+mt for mt in mutTypes]
    sampleTimepointDF = pd.DataFrame(sampleTimepointDF, columns=sampleTimepointDFcolumns)

    ratesDFrowList = [] # list to be populated with correlations between all types of mutations and timepoints, yielding a mutation rate as slope
    ratesDFrowColumns = ['sampleLabel', 'replicate', 'rate']

    # loop through each row in the timepoints table, then loop through each column in the timepoints table,
    #   grab the total mutations as well as the total of each mutation type, divide by # of sequences,
    #   background subtract if background barcode group is given in config file, divide by # of bases per sequence
    #   (dependent on both sequence length and identity), yielding mutations per base analyzed of each type

    # loop through each row in timepoints, then loop through each column in the timepoints table
    for sampleLabel, row in timepointsCSV.iterrows():

        for generations in timepointsCSV.columns:
            timepointTag, timepointBCgroup = row[generations].split('_')
            timepointRefSeqfasta = config['runs'][timepointTag]['reference']
            timepointRefSeq = str(list(SeqIO.parse(timepointRefSeqfasta, 'fasta'))[1].seq).upper()

            # grab row from mut stats corresponding to sample/barcode group
            timepointSeqsMutStatsRow = (mutStatsCSV.loc[(mutStatsCSV['tag']==timepointTag) & (mutStatsCSV['barcode_group']==timepointBCgroup)]).iloc[0]

            # calculate the substitutions per base analyzed for all types of substitutions, normalized to # of sequences
            normTrimmedRow = trim_normalize_row(timepointSeqsMutStatsRow, timepointRefSeq, mutTypes)
            if backgroundBool:
                normTrimmedRow = list(np.array(normTrimmedRow) - np.array(backgroundRows[timepointTag]))
            sampleTimepointDFrowList.append([sampleLabel, generations] + normTrimmedRow)
        
    sampleTimepointDF = pd.DataFrame(sampleTimepointDFrowList, columns=sampleTimepointDFcolumns)


    print(sampleTimepointDF)

    plotList = []

    # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)

    # generate random data
    df = pd.DataFrame(data={'col_1': np.random.normal(5, 2, 100)})

    df['col_2'] = df.col_1 + np.random.gamma(5, 2, 100)
    df['col_3'] = df.col_1*2 + np.random.normal(0, 10, 100)
    df['col_4'] = df.col_1**2 + np.random.normal(0, 10, 100)
    df['col_5'] = np.sin(df.col_1)
    df['col_6'] = np.cos(df.col_1)
    corr = df.corr().abs()
    # mask the upper triangle of the heatmap
    corr.values[np.triu_indices_from(corr, 0)] = np.nan

    heatmap = hv.HeatMap((corr.columns, corr.index, corr))\
                .opts(tools=['hover'],  height=400, width=400, fontsize=9,
                    toolbar='above', colorbar=False, cmap='Blues',
                    invert_yaxis=True, xrotation=90, xlabel='', ylabel='',
                    title='Correlation Coefficient Heatmap (absolute value)')

    # define tap stream with heatmap as source
    tap_xy = hv.streams.Tap(source=heatmap, x='col_1', y='col_4')

    # calculate correlation plot based on tap
    def tap_corrplot(x, y):
        # drop missing values if there are any
        df_notnull = df[[x, y]].dropna(how='any')

        # fit a 2nd degree line/curve
        m1, m2, b = np.polyfit(df_notnull[x], df_notnull[y], deg=2)
        # generate data to plot fitted line/curve
        x_curve = np.linspace(df[x].min(), df[x].max())
        y_curve = m1*x_curve**2 + m2*x_curve+ b

        curve = hv.Curve((x_curve, y_curve), x, y)\
                .opts(color='#fc4f30', framewise=True)

        scatter = hv.Scatter((df[x], df[y]), x, y)\
                    .opts(height=400, width=400, fontsize=9, size=5,
                        alpha=0.2, ylim=(df[y].min(), df[y].max()),
                        color='#30a2da', framewise=True,
                        title='Correlation Plot (2nd degree fit)')

        return curve * scatter

    # dynamic map cant be saved as .html, so cant use
    # maybe can use pd.DataFrame.corr()

    # map tap in heatmap with correlation plot
    tap_dmap = hv.DynamicMap(tap_corrplot, streams=[tap_xy])

    p = heatmap + tap_dmap

    # defaults = dict(width=800, height=800, xaxis=None, yaxis=None, tools=['tap', 'hover', 'box_select'])

    # p.opts(
    #     hv.opts.Scatter( color_index=2, tools=['tap', 'hover'], width=600, framewise=True, marker='triangle', cmap='Set1', size=10),
    #     hv.opts.Overlay( toolbar='above', legend_position='right' ),
    #     hv.opts.Curve( line_color='black', framewise=True ) )

    # def gen_samples(N, corr=0.8):
    #     xx = np.array([-0.51, 51.2])
    #     yy = np.array([0.33, 51.6])
    #     means = [xx.mean(), yy.mean()]  
    #     stds = [xx.std() / 3, yy.std() / 3]
    #     covs = [[stds[0]**2          , stds[0]*stds[1]*corr], 
    #             [stds[0]*stds[1]*corr,           stds[1]**2]] 

    #     return np.random.multivariate_normal(means, covs, N)

    # data = [('Week %d' % (i%10), np.random.rand(), chr(65+np.random.randint(5)), i) for i in range(100)]
    # sample_data = hv.NdOverlay({i: hv.Points(gen_samples(np.random.randint(1000, 5000), r2))
    #                             for _, r2, _, i in data})
    # points = hv.Scatter(data, kdims=['Date', 'r2'], vdims=['block', 'id']).redim.range(r2=(0., 1))
    # stream = Selection1D(source=points)
    # empty = (hv.Points(np.random.rand(0, 2)) * hv.Curve(np.random.rand(0, 2))).relabel('No selection')

    # def regression(index):
    #     if not index:
    #         return empty
    #     scatter = sample_data[index[0]]
    #     xs, ys = scatter['x'], scatter['y']
    #     slope, intercep, rval, pval, std = stats.linregress(xs, ys)
    #     xs = np.linspace(*scatter.range(0)+(2,))
    #     reg = slope*xs+intercep
    #     return (scatter * hv.Curve((xs, reg))).relabel('r2: %.3f' % slope)

    # reg = hv.DynamicMap(regression, kdims=[], streams=[stream])

    # average = hv.Curve(points, kdims=['Date'], vdims=['r2']).aggregate(function=np.mean)
    # p = points * average + reg

    hv.save(p, snakemake.output.rate, backend='bokeh')

    # raise(RuntimeError)



    # if True:

    #     timepointSamples[generations] = tuple(row[generations].split('_'))

    #     if row['total_seqs'] == 0: continue

    #     bcGroupSpectrum = normalized_spectrum_df(row, backgroundRow, backgroundBool, mutTypes, nt_normal_dict)

    #     bcGroupSpectrumMelted = pd.melt(bcGroupSpectrum, value_vars=bcGroupSpectrum)
    #     bcGroupSpectrumMelted['wtNT'] = bcGroupSpectrumMelted.apply(lambda row:
    #         row['variable'][0], axis=1)
    #     bcGroupSpectrumMelted['mutNT'] = bcGroupSpectrumMelted.apply(lambda row:
    #         row['variable'][3], axis=1)
        
    #     bcGroupSpectrumPivot = bcGroupSpectrumMelted.pivot(index='wtNT', columns='mutNT', values = 'value')

    #     dataOutName = os.path.join(csvOutDir, f'{tag}_{row.barcode_group}_mutSpectrum.csv')
    #     renameDict = {}
    #     for nt in list('ATGC'): renameDict[nt] = 'mut: '+nt
    #     bcGroupSpectrumPivot.rename(columns=renameDict).to_csv(dataOutName)

    #     bcGroupSpectrumPivot = bcGroupSpectrumPivot.fillna(0.0) # replace NaNs with 0
    #     bcGroupSpectrumPivot = bcGroupSpectrumPivot.clip(lower=0) # convert negative values to 0
        
    #     TOOLTIPS = [('mutation type', '@wtNT'+'->'+'$name'), ('proportion of mutations', '@$name')]
    #     plotTitle = f"{tag}_{row['barcode_group']}"
    #     spectrumPlot = figure(title=plotTitle, plot_height=400, plot_width=400, x_range=list('ATGC'), tooltips=TOOLTIPS)

    #     if first: # legend only for first plot
    #         spectrumPlot.vbar_stack(list('ATGC'), x='wtNT', width = 0.9, source=ColumnDataSource(bcGroupSpectrumPivot),
    #             color = Blues5[:4], legend_label=list('ATGC'))
    #         spectrumPlot.legend.title = 'mutation'
    #         first = False
    #     else:
    #         spectrumPlot.vbar_stack(list('ATGC'), x='wtNT', width = 0.9, source=ColumnDataSource(bcGroupSpectrumPivot),
    #             color = Blues5[:4])
        
    #     spectrumPlot.yaxis.axis_label = 'normalized substitution frequency'

    #     plotList.append(spectrumPlot)
    # if len(plotList)>0:
    #     plotList[-1].xaxis.axis_label = 'wild type nucleotide' # label x axis for bottom plot

    # save(column(plotList))

if __name__ == '__main__':
    main()