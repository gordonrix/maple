#
#  DESCRIPTION   : Script for maple pipeline. Generates a scatter plot of genotypes with
#                       point size and color determined by user-defined columns
#
#  AUTHOR(S)     : Gordon Rix
#

import hvplot.pandas
import pandas as pd
import colorcet as cc
from pandas.api.types import is_numeric_dtype

data = pd.read_csv(snakemake.input.genotypesReduced)
data = data[data.loc[:,['dim1','dim2']].notnull().all(axis=1)] # ignore genotypes that could not be assigned x/y values because of indels
size_column = snakemake.params.size_column
color_column = snakemake.params.color_column

if not all([x in list(data.columns) for x in [size_column, color_column]]):
    print("[WARNING] Only columns found in the genotypes.csv file may be used to map point size and color")

minSize, maxSize = [int(x) for x in snakemake.params.size_range.replace(' ','').split(',')]
if not minSize<=maxSize:
    (f"For genotypes2D plot, minimum size must be less than maximum size min/max of {minSize}/{maxSize} provided")

# assign point size
maxSizeCol = data[size_column].max()
minSizeCol = data[size_column].min()
if maxSizeCol == minSizeCol:
    data.loc[:,'point_size'] = minSize # use mininum point size for all
else:
    data['point_size'] = data[size_column].apply( lambda x:
        ((x-minSizeCol) / (maxSizeCol-minSizeCol)) * (maxSize - minSize) + minSize)

# assign color
colorDict = {} # dictionary of value:color key:value pairs to map color_column variables to specific colors
if is_numeric_dtype(data[color_column]):                            # color by value for numerical column
    legendBool = True
    colormap = cc.blues
    minColorCol, maxColorCol = data[color_column].min(), data[color_column].max()
    colorConstant = (len(colormap)-1) / ( maxColorCol - minColorCol )
    for ccValue in data[color_column].unique():
        colorDict[ccValue] = colormap[int((ccValue-minColorCol)*colorConstant)]
else:                                                               # random colors for non-numerical column
    legendBool = False
    colormap = cc.bmy   
    colorConstant = len(colormap) / len(data[color_column].unique())
    for i, ccValue in enumerate(data[color_column].unique()):
        colorDict[ccValue] = colormap[int(i*colorConstant)]
data['color'] = data[color_column].map(colorDict)

plot = data.hvplot.scatter(x='dim1', y='dim2', size='point_size', by=color_column, color='color', legend=legendBool, hover_cols=list(data.columns)[:10], width=1200, height=1000).opts(
    xaxis=None, yaxis=None)
hvplot.save(plot, snakemake.output.genotypes2Dplot)