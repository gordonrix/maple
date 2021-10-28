from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import column
from bokeh.models import ColumnDataSource
import os
import pandas as pd
import numpy as np
from snakemake.io import Namedlist

### Asign variables from config file and inputs
config = snakemake.config
tag = snakemake.wildcards.tag
AAorNT = snakemake.wildcards.AAorNT
N = config[f'{AAorNT}_distribution_plot_x_max']
###

### Output variables
output_file(snakemake.output[0])
###

if type(snakemake.input.dist) == Namedlist:
    mode = 'grouped'
    inputList = snakemake.input.dist
else:
    mode = 'individual'
    inputList = [snakemake.input.dist]

inputList = sorted(inputList)

plotDict = {}

for inFile in inputList:    

    bc = inFile.split('_')[-2]
    data = pd.read_csv(inFile)
    # if int(data.iloc[:,[1]].sum())==0: continue
    
    plotTitle = f"{tag}_{bc}"

    xLabel, yLabel = data.columns

    TOOLTIPS = [(xLabel, f'@{xLabel}'), (yLabel, f'@{yLabel}')]
    plotDict[bc] = figure(title=plotTitle, plot_width=600, plot_height=400, x_range=(-0.7,N), tooltips=TOOLTIPS)
    plotDict[bc].vbar(x=xLabel, top=yLabel, width=0.5, bottom=0, color='black', source=ColumnDataSource(data))
    plotDict[bc].xaxis.axis_label = xLabel
    plotDict[bc].yaxis.axis_label = yLabel

save(column([plotDict[key] for key in plotDict]))
