from bokeh.plotting import figure, output_file, show, save
from bokeh.layouts import column
import os
import pandas as pd
import numpy as np
from snakemake.io import Namedlist

### Asign variables from config file and inputs
config = snakemake.config
tag = snakemake.wildcards.tag
AAorNT = snakemake.wildcards.AAorNT
N = config['percentile']
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

plotDict = {}

def calculate_Nth_percentile(N, DF):
    """
    Given distribution as a dataframe (DF), calculates (N)th percentile
    of the distribution
        ex: calculate_Nth_percentile(N, dataframe) returns value (maxX) at which
        90% of distribution is less than maxX
    """
    distSeries = DF[f'seqs_with_n_{AAorNT}substitutions']
    numSeqsInNth = (int(round((N/100)*distSeries.sum())))
    runningTotalSeqs = 0
    for i, numSeqs in enumerate(distSeries):
        runningTotalSeqs += numSeqs
        maxX = i
        if runningTotalSeqs > numSeqsInNth:
            break
    return maxX

distInfoDictList = []
percentiles = []
for inFile in inputList:    # get the maximum Nth percentile for all plots

    distInfoDict = {}

    distInfoDict['barcodes'] = inFile.split('_')[-2]
    distInfoDict['df'] = pd.read_csv(inFile)

    if int(distInfoDict['df'].iloc[:,[1]].sum())==0: continue

    distInfoDictList.append(distInfoDict)

    percentiles.append(calculate_Nth_percentile(N, distInfoDict['df']))

maxpercentile = np.max(percentiles)
if maxpercentile < 5:
    maxpercentile = 5

for distInfoDict in distInfoDictList:
    
    plotTitle = f"{tag}_{distInfoDict['barcodes']}"

    # lastNonZero = (dist.shape[0]-dist.ne(0).values[::-1].argmax(0)-1)[1]
    distTrimmed = distInfoDict['df'][:int(maxpercentile*1.5)]
    xLabel, yLabel = distTrimmed.columns

    bc = distInfoDict['barcodes']
    plotDict[bc] = figure(title=plotTitle, plot_width=600, plot_height=400)
    plotDict[bc].vbar(x=distTrimmed['n'], top=distTrimmed.iloc[:,1], width=0.5, bottom=0, color='black')
    plotDict[bc].xaxis.axis_label = xLabel
    plotDict[bc].yaxis.axis_label = yLabel

save(column([plotDict[key] for key in plotDict]))
