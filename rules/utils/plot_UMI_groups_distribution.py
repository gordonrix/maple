"""part of nanopype-MACE pipeline, written by Gordon Rix
plot_UMI_groups_distribution.py
plots the distribution of counts of UMIs identified in provided UMI_tools group log file"""

from bokeh.plotting import figure, output_file, show, save
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

# get data, convert into a proper distribution
fullDF = pd.read_csv(inFile, sep='\t')
UMIcounts = fullDF[['final_umi_count', 'unique_id']].drop_duplicates()['final_umi_count']
UMIcountsList = list(UMIcounts)
maxCount = np.max(UMIcountsList)
dist = np.zeros([maxCount])
for count in UMIcountsList:
    dist[count-1] += 1

xLabel, yLabel = 'n', 'unique UMI groups with n sequences'
plotDF = pd.DataFrame({xLabel:[i for i in range(1,maxCount+1)], yLabel:list(dist)})

plot = figure(title=plotTitle, plot_width=800, plot_height=600)
plot.vbar(x=plotDF[xLabel], top=plotDF.iloc[:,1], width=0.5, bottom=0, color='black')
plot.xaxis.axis_label = xLabel
plot.yaxis.axis_label = yLabel

plotDF.to_csv(outCSV)
save(plot)