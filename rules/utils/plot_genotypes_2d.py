#
#  DESCRIPTION   : Script for maple pipeline. Generates a scatter plot of genotypes with
#                       point size and color determined by user-defined columns
#
#  AUTHOR(S)     : Gordon Rix
#

import hvplot.pandas
import holoviews as hv
import pandas as pd
import numpy as np
import colorcet as cc
from pandas.api.types import is_numeric_dtype
from common import export_svg_plots, cmap_dict

data = pd.read_csv(snakemake.input.genotypesReduced)

if snakemake.params.plot_AA:
    dim1,dim2 = 'AA_PaCMAP1', 'AA_PaCMAP2'
else:
    dim1,dim2 = 'NT_PaCMAP1', 'NT_PaCMAP2'

data = data[data.loc[:,[dim1,dim2]].notnull().all(axis=1)] # ignore genotypes that could not be assigned x/y values because of indels

downsample = snakemake.params.downsample
if downsample:
    if downsample < len(data):
        idx = np.sort(np.random.choice(np.arange(len(data)), size=downsample, replace=False))
        data = data.iloc[idx]

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
    legendBool = False
else:                                                               # random colors for non-numerical column
    legendBool = True

colormap = cmap_dict()[snakemake.params.cmap]

scatter = data.hvplot(kind='points', x=dim1, y=dim2, size='point_size', color=color_column, clabel=color_column,
    cmap=colormap, legend=legendBool, xticks=[100], yticks=[100]).opts(
    tools=[], xlabel=dim1, ylabel=dim2, height=800, width=905)  # slightly wider to account for colorbar so that plot dimensions match hexbins

hexbins = data.hvplot.hexbin(x=dim1, y=dim2, clabel='Count', color=color_column,
    cmap='dimgray', xticks=[100], yticks=[100]).opts(
    xlabel=dim1, ylabel=dim2, height=800, width=800)

try:
    export_svg_plots([scatter, hexbins], snakemake.output.genotypes2Dscatter, labels=['scatter','hexbins'])
except:
    print(f"""[ERROR] SVG export for {snakemake.output.genotypes2Dscatter} failed. Pipeline continuing but the SVG version of this plot was not generated. 
            Try again by deleting the associated .html file and rerunning snakemake.""")

hvplot.save(scatter, snakemake.output.genotypes2Dscatter)
hvplot.save(hexbins, snakemake.output.genotypes2Dhexbins)

