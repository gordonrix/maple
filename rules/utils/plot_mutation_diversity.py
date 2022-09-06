""" script for maple pipeline

Visualization of hamming distances as a distribution
bar plot and a network graph of unique genotypes
"""

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import column
from bokeh.models import ColumnDataSource
import os
import math
import pandas as pd
import numpy as np
import networkx as nx
import holoviews as hv
from bokeh.io import output_file, save
from yaml.nodes import SequenceNode
from Bio import SeqIO
from Bio.Seq import translate

### Asign variables from config file
config = snakemake.config
tag = snakemake.wildcards.tag
bc = snakemake.wildcards.barcodes
maxHD = config.get('diversity_plot_hamming_distance_edge_limit', False)
###

def plot_distribution(tag, bc, title, binCountsDF, outputName, xmax):
    """
    generate and export a distribution plot from a dataframe of bincounts"
    """
    plotTitle = f"{title}, {tag}_{bc}"
    xLabel, yLabel = binCountsDF.columns
    TOOLTIPS = [(xLabel, f'@{xLabel}'), (yLabel, f'@{yLabel}')]
    hamDistPlot = figure(title=plotTitle, plot_width=600, plot_height=400, x_range=(-0.7, xmax), tooltips=TOOLTIPS)
    hamDistPlot.vbar(x=xLabel, top=yLabel, width=0.5, bottom=0, color='black', source=ColumnDataSource(binCountsDF))
    hamDistPlot.xaxis.axis_label = xLabel
    hamDistPlot.yaxis.axis_label = yLabel
    output_file(outputName)
    save(hamDistPlot)

if len(snakemake.input) == 1:
    aaHDdist = pd.read_csv(snakemake.input.aaHamDistCSV)
    plot_distribution(tag, bc, "AA pairwise hamming distances", aaHDdist, snakemake.output.aaHamDistPlot, config['hamming_distance_distribution_plot_x_max'])

else: # make nucleotide hamming distance and network graph

    genotypesDF = pd.read_csv(snakemake.input.genotypes, dtype={'genotype_ID':str}, na_filter=False)
    genotypesDF.drop(['NT_insertions','NT_deletions'], axis=1, inplace=True)
    hammingDistanceEdgesDF = pd.read_csv(snakemake.input.edges, dtype={'source':str,'target':str})
    # only use genotypes present in the edges DF, which only includes genotypes that have been filtered by mutation_diversity.py
    filteredGenotypes = hammingDistanceEdgesDF['source'].unique().tolist() + hammingDistanceEdgesDF['target'].unique().tolist()
    genotypesDF = genotypesDF[genotypesDF['genotype_ID'].isin(filteredGenotypes)]
    genotypesDF.set_index('genotype_ID', inplace=True)

    ntHDdist = pd.read_csv(snakemake.input.ntHamDistCSV)
    plot_distribution(tag, bc, "NT pairwise hamming distances", ntHDdist, snakemake.output.ntHamDistPlot, config['hamming_distance_distribution_plot_x_max'])

    ### generate network graph of sequences as nodes and edges connecting all nodes, with inverse hamming distance as edge weight. Plot with holoviews
    if snakemake.params.edgeLimit:
        hammingDistanceEdgesDF = hammingDistanceEdgesDF.loc[hammingDistanceEdgesDF['hammingDistance']<=snakemake.params.edgeLimit]                                       # apply filter for max hamming distance based on user-supplied value
    else: 
        hammingDistanceEdgesDF = hammingDistanceEdgesDF[hammingDistanceEdgesDF['hammingDistance'] < max(3, np.argmax(ntHDdist))]    # filter out edges with hamming distance greater than or equal to (a) the maximum hamming distance bincount (will be median for normal distribution), or (b) 3, whichever is larger
    if len(hammingDistanceEdgesDF)==0:
        exit('[ERROR] plot_mutation_diversity failed because there were no sequence pairs to process after applying the hamming distance cutoff.')
    def mutCountHDweighting(source,target, hammingDistance):
        sourceMutCount = int(genotypesDF.at[source,'NT_substitutions_count'])
        targetMutCount = int(genotypesDF.at[target,'NT_substitutions_count'])
        return np.e**((sourceMutCount+targetMutCount)/2-hammingDistance)       # weight is e^(the average mutation count of the two mutants minus the hamming distance). This makes weight dependent on both hamming distance and # of mutations, resulting in better clustering of similar clades

    hammingDistanceEdgesDF['weight'] = hammingDistanceEdgesDF.apply(lambda row:
        mutCountHDweighting(row['source'], row['target'], row['hammingDistance']), axis=1)
    # hammingDistanceEdgesDF['weight'] = 2**(3-hammingDistanceEdgesDF['hammingDistance'])       # weighting based on hamming distance alone. Result is less structured plot that just has concentric rings of nodes that track with increasing hamming distance from WT
    G = nx.from_pandas_edgelist(hammingDistanceEdgesDF, edge_attr='weight')
    nx.set_node_attributes(G, genotypesDF.to_dict(orient='index'))
    nx.write_gexf(G, snakemake.output.GraphFile)
    hv.extension('bokeh')
    defaults = dict(width=800, height=800, xaxis=None, yaxis=None, tools=['tap', 'hover', 'box_select'])
    hv.opts.defaults(
        hv.opts.EdgePaths(**defaults), hv.opts.Graph(**defaults), hv.opts.Nodes(**defaults))
    nodeColorMap = 'blues' if config['force_directed_plot_node_color'] in ['count', 'NT_substitutions_count', 'AA_substitutions_nonsynonymous_count', ] else 'rainbow'

    # linear equation to scale Node sizes to range from 5-25 (for node weight)
    maxCount = genotypesDF[config['force_directed_plot_node_size']].max()
    minCount = genotypesDF[config['force_directed_plot_node_size']].min()
    if maxCount==minCount: # set all to 10 if no difference between max and min
        slopeN, interceptN = 0, 10
    else:
        slopeN = (25-5) / (maxCount-minCount)
        interceptN = 5 - (slopeN*minCount)

    # equation to scale log of edge widths 0.05-7 (for edge weight)
    maxWeight = np.log(hammingDistanceEdgesDF['weight'].max())
    minWeight = np.log(hammingDistanceEdgesDF['weight'].min())
    if maxWeight==minWeight: # set all to 1 if no difference between max and min
        slopeW, interceptW = 0, 1
    else:
        slopeW = (10-0.05) / (maxWeight-minWeight)
        interceptW = 0.1 - (slopeW*minWeight)

    # equation to scale log of edge widths 0.05-0.4 (for edge Alpha, or opacity)
    if maxWeight==minWeight: # set all to 0.4 if no difference between max and min
        slopeA, interceptA = 0, 0.4
    else:
        slopeA = (0.4-0.05) / np.absolute(maxWeight-minWeight)
        interceptA = 0.1 - (slopeA*minWeight)

    networkPlot = hv.Graph.from_networkx(G, nx.layout.fruchterman_reingold_layout).opts(
        hv.opts.Graph(node_size=(hv.dim(config['force_directed_plot_node_size'])*slopeN)+interceptN, node_color=config['force_directed_plot_node_color'], cmap=nodeColorMap,
                        edge_line_width=(np.log(hv.dim('weight'))*slopeW)+interceptW, edge_color=np.log(hv.dim('weight')), edge_cmap='Inferno', edge_alpha=np.log(hv.dim('weight'))*slopeA+interceptA))
    hv.save(networkPlot, snakemake.output.GraphPlot, backend='bokeh')