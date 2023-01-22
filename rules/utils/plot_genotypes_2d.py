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
from bokeh.models import HoverTool

data = pd.read_csv(snakemake.input.genotypesReduced)
data = data[data.loc[:,['dim1','dim2']].notnull().all(axis=1)] # ignore genotypes that could not be assigned x/y values because of indels

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
    colormap = 'kbc_r'
else:                                                               # random colors for non-numerical column
    legendBool = False
    colormap = 'bmy'

hover = HoverTool(tooltips=[('count','@count'),('NT mutations count','@NT_substitutions_count'),('AA mutations count','@AA_substitutions_nonsynonymous_count'),
                            ('NT mutations','@NT_substitutions'),('AA mutations','@AA_substitutions_nonsynonymous')])
tools = ['box_select', 'lasso_select',hover]
plot = data.hvplot(kind='points', x='dim1', y='dim2', size='point_size', color=color_column, hover_cols=[color_column, 'count', 'NT_substitutions_count', 'AA_substitutions_nonsynonymous_count', 'NT_substitutions', 'AA_substitutions_nonsynonymous'],
    legend=legendBool, width=1000, height=800).opts(
    xaxis=None, yaxis=None)

hvplot.save(plot, snakemake.output.genotypes2Dplot)