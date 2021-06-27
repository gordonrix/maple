import os
import random
from math import pi

import numpy as np
import pandas as pd
from bokeh.layouts import column
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource, FactorRange,
                          HoverTool, Legend, LinearColorMapper,
                          PrintfTickFormatter)
from bokeh.palettes import Inferno
from bokeh.plotting import figure, output_file, save, show
from snakemake.io import Namedlist

### Asign variables from config file and inputs
config = snakemake.config
tag = snakemake.wildcards.tag
AAorNT = snakemake.wildcards.AAorNT
inputFrequencies = snakemake.input.frequencies
mutStatsDF = pd.read_csv(snakemake.input.mutStats, dtype={'tag':str,'barcode_group':str})
if config['mutations_frequencies_raw'] == True:
    yAxisLabel = 'total mutation count'
elif config['mutations_frequencies_raw'] == False:
    yAxisLabel = 'proportion of reads'
###

output_file(snakemake.output[0])

if type(inputFrequencies) == Namedlist:
    mode = 'grouped'
else:
    mode = 'individual'
    inputFrequencies = [inputFrequencies]

wtColumn = f'wt_{AAorNT}'

plotList = []
first = True #add legend only for first plot
for inFile in inputFrequencies:

    ### Get mutation data, convert to more readily plottable format
    barcodeGroup = inFile.split('_')[-2]
    plotTitle = f'{tag}_{barcodeGroup}'
    mutsDF = pd.read_csv(inFile, index_col=0).transpose()
    wt = list(mutsDF.index)
    totalSequences = mutStatsDF.loc[mutStatsDF['barcode_group']==barcodeGroup, 'total_seqs'].iloc[0]

    if totalSequences == 0:
        continue

    # set appropriate dataframe and tooltips label to use for plotting
    if yAxisLabel == 'total mutation count':
        tooltipValueLabel = 'total count'
    elif yAxisLabel == 'proportion of reads':
        tooltipValueLabel = 'proportion of reads'

    yAxisLabelWithReadCount = yAxisLabel + f' (n = {totalSequences})'

    # setting colors and hatch patterns
    if AAorNT == 'AA':

        palette = Inferno[6]
        newPalette = []
        for p in palette[::-1]:
            newPalette.append(p)
        palette = newPalette

        muts =     ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W',
                    'S', 'T', 'N', 'Q',
                    'R', 'H', 'K',
                    'D', 'E',
                    'C', 'G', 'P',
                    '*']

        mutsDF = mutsDF[muts]

        colors =  []    # assign colors to different mutation types
        for i,n in enumerate([8, 4, 3, 2, 3, 1]):
            colors.extend([palette[i] for x in range(0,n)])

        hatches =  [' ', 'v', '"', '>', '.', '`', 'x', 'o',
                    ' ', 'v', '"', '>',
                    ' ', 'v', '"',
                    ' ', 'v',
                    ' ', 'v', '"',
                    ' ']

        hatchColors = []
        for c,n in zip(['#000000', '#000000', '#D4DCDE', '#D4DCDE', '#D4DCDE', '#D4DCDE'],[8, 4, 3, 2, 3, 1]):
            hatchColors.extend([c for x in range(0,n)])

    elif AAorNT == 'NT':
        colors = Inferno[4]
        hatches = [' ', ' ' ,' ' ,' ']
        hatchColors = ['#000000', '#000000', '#000000', '#000000']
        muts = list(mutsDF.columns)

    pHeight = 250
    if first:
        pHeight = int(pHeight*1.3)
        maxCount = max(list(mutsDF.max(axis=0)))

    mutsDF = mutsDF.reset_index().rename(columns={'index':wtColumn})
    source = ColumnDataSource(mutsDF)

    TOOLTIPS = [('mutation', f'@{wtColumn}'+'$name'), (tooltipValueLabel, '@$name')]
    barcodePlot = figure(title=plotTitle, plot_height=pHeight, plot_width=len(mutsDF)*12, x_range=FactorRange(*wt), tooltips=TOOLTIPS)

    v = barcodePlot.vbar_stack(muts, x=wtColumn, width=0.9, color=colors, source=source,
        legend_label=muts, hatch_pattern=hatches, hatch_scale=5.0, hatch_color=hatchColors)
    barcodePlot.xaxis.major_label_orientation = pi/2.6
    barcodePlot.xaxis.axis_label = f'wild type {AAorNT}'
    barcodePlot.yaxis.axis_label = yAxisLabelWithReadCount

    # add legend only for first (top) plot
    if first:
        barcodePlot.y_range.end = maxCount*1.5 # make first plot taller to accomodate legend
        barcodePlot.legend.location = 'top_left'
        barcodePlot.legend.orientation = 'horizontal'
        first = False
    else:
        barcodePlot.legend.visible = False

    plotList.append(barcodePlot)

save(column(plotList))
