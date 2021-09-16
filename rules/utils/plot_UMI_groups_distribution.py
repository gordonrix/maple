"""part of nanopype-MACE pipeline, written by Gordon Rix
plot_UMI_groups_distribution.py
plots the distribution of counts of UMIs identified in provided UMI_tools group log file"""

from bokeh.plotting import figure, output_file, show, save
from bokeh.models import (BasicTicker, ColorBar, ColumnDataSource, FactorRange,
                          HoverTool, Legend, LinearColorMapper,
                          PrintfTickFormatter)
import os
import pandas as pd
import numpy as np

### Asign variables from config file, inputs, and outputs
config = snakemake.config
tag = snakemake.wildcards.tag
inFile = snakemake.input[0]
outCSV = snakemake.output.csv
output_file(snakemake.output.plot)
###
    
plotTitle = f"{tag}_UMI_groups_distribution"

# get data, convert into a proper distribution, with one column representing number of read groups, and one representing total reads
fullDF = pd.read_csv(inFile, sep='\t')
UMIcounts = fullDF[['final_umi_count', 'unique_id']].drop_duplicates()['final_umi_count']
UMIcountsList = list(UMIcounts)
maxCount = np.max(UMIcountsList)
dist = np.zeros([maxCount])
for count in UMIcountsList:
    dist[count-1] += 1

xLabel, yLabel = 'n_number_of_reads_in_UMI_group', 'reads_in_UMI_groups_with_n_reads'
# xLabel, yLabel = 'x', 'y'
plotDF = pd.DataFrame({xLabel:[i for i in range(1,maxCount+1)], yLabel:[(i+1)*a for i,a in enumerate(dist)], 'unique_UMI_groups':list(dist)})

TOOLTIPS = [(xLabel, f'@{xLabel}'), (yLabel, f'@{yLabel}'), ('unique_UMI_groups', '@unique_UMI_groups')]
# TOOLTIPS = [(xLabel, '$name'), (yLabel, '$name'), ('unique UMI groups', '$name')] 
plot = figure(title=plotTitle, plot_width=800, plot_height=600, tooltips=TOOLTIPS)
plot.vbar(x=xLabel, top=yLabel, width=0.5, bottom=0, color='black', source=ColumnDataSource(plotDF))
plot.xaxis.axis_label = xLabel
plot.yaxis.axis_label = yLabel

plotDF.to_csv(outCSV, index=False)
save(plot)