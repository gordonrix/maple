""" script for maple pipeline

Uses output data files from rule mutation_analysis for all files being processed,
calculates statistics on diversity and visualizes this information as a distribution
of hamming distances and a network graph of sequences
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

### Asign variables from config file
config = snakemake.config
tag = snakemake.wildcards.tag
bc = snakemake.wildcards.barcodes
###

genotypesDF = pd.read_csv(snakemake.input[0], na_filter=False)
genotypesDF.drop(genotypesDF[genotypesDF['NT_insertions'] != ''].index, inplace=True)
genotypesDF.drop(genotypesDF[genotypesDF['NT_deletions'] != ''].index, inplace=True)
genotypesDF.drop(genotypesDF[genotypesDF['count'] == 0].index, inplace=True) # removes wildtype row if count is 0
genotypesDF.drop(['NT_insertions','NT_deletions'], axis=1, inplace=True)
genotypesDF.set_index('genotype', inplace=True)
refSeqfasta = config['runs'][tag]['reference']
NTref = list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq
NTrefLength = len(NTref)
NTs = 'ATGC'

# construct a 3 dimensional array of all sequences encoded as 2 dimensional arrays, which are the length of the sequence and all four possible nucleotides, and the third dimension
#   is the genotype ID of the sequence in the order given by the genotypes .csv file, with wildtype (an array of 0s) as the first sequence. From this 3d array,
#   a 2d array will be contructed, consisting of the pairwise hamming distances between each pair of sequences
seqLength = len(NTref)
numNTs = len(NTs)
numGenotypes = len(genotypesDF)
genotypes3DArray = np.zeros((seqLength, numNTs, numGenotypes), dtype=int)

def sequence_array_from_NT_substitutions(NTsubstitutions, zeroArray):
    """
    from a str-type list of nucleotide substitutions and a numpy zero array of shape = (length of nucleotide sequence, 4),
    will flip any 0s corresponding to a mutation to 1 and return the resulting array
    """
    if NTsubstitutions == '':
        return zeroArray
    NTsubstitutionsList = NTsubstitutions.split(', ')
    count = 0
    for mutation in NTsubstitutionsList:
        mutNT = mutation[-1]
        posi = int(mutation[1:-1])
        zeroArray[posi,NTs.find(mutNT)] += 1
        count += 1
    return zeroArray

for i, subsList in enumerate(genotypesDF['NT_substitutions']):
    zeroArray = np.zeros((int(len(NTref)), len(NTs)), dtype=int)
    genotypes3DArray[:,:,i] = sequence_array_from_NT_substitutions(subsList, zeroArray)

hammingDistanceDF = np.empty((len(genotypesDF), len(genotypesDF)), dtype=int) * np.nan

matrixPad = 0 # number of NaNs to pad the beginning of each row so that the hamming distance positions match up with the columns. These are essentially placeholders for hamming distances that were already calculated, and can be found on the other half of the N by N square
hammingDistanceRows = []
hammingDistanceBinCounts = np.zeros(NTrefLength, dtype=int)    # Hamming distances for each row converted into bin counts and multiplied by # of times sequence is observed in the dataset

# iterate through each row and generate counts of hamming distances (bincounts) based on hamming distance between the row genotype and
#   all other genotypes, including itself, to be used for plotting distribution, and also generate a matrix of hamming distances between
#   the row genotype and all other genotypes, excluding itself, to be used for plotting a network graph
for row in range(0, numGenotypes):
    matrixPad += 1

    hammingDistance3Darray = np.absolute( genotypes3DArray[:,:,row:] - genotypes3DArray[:,:,row].reshape(seqLength, numNTs, 1) )    # subtract the 3D array representing remaining genotypes from 3D array representing many copies of the genotype for this row, then take absolute value of all elements in this 3D array
    HDrow = np.sum(hammingDistance3Darray, axis=1)                                                                                  # sum across nucleotides to get 2D array of total mutations at each position and for each genotype
    HDrow[HDrow>1] = 1                                                                                                              # set values above 1 to 1. otherwise, genotypes that both have different mutations at a single position will result in a hamming distance of 2 for that position, but hamming distance can only ever be 1 for a single position
    HDrow = np.sum(HDrow, axis=0)                                                                                                   # sum across position axis to give 1D array of total hamming distance for each genotype
    hammingDistanceRows.append( np.pad(HDrow[1:], ((matrixPad,0)), constant_values=-1) )                                            # make new row of all hamming distances for the row genotype, padding the beginning of the row with the null value -1 so that when combined the column positions will line up

    genotypeCount = int(genotypesDF.iloc[[row]]['count'])
    HDcountsList = [0] * int((genotypeCount*(genotypeCount-1))/2)                                                          # 0 count for comparison of the row genotype to itself calculated as n(n-1)/2, see https://stackoverflow.com/questions/18859430/how-do-i-get-the-total-number-of-unique-pairs-of-a-set-in-the-database

    # loop through hamming distances for row except first one which is always 0
    for HDindex in range(1, len(HDrow)):

        # increase counts of each hamming distance by the product of the counts for the two sequences being compared
        HDcountsList.extend( [HDrow[HDindex]] * int(genotypesDF.iloc[[row]]['count']) * int(genotypesDF.iloc[[HDindex+row]]['count']) )

    HDrowBins = np.bincount(HDcountsList)                                                                               # get counts of each hamming distance observed, then multiply these bincounts by the count for the row genotype to give the total count of these observed hamming distances
    HDrowBins.resize(NTrefLength)                                                                                       # reshape to be the same length as the sequence, which is the maximum possible hamming distance, allowing for concatenation with bincounts for other genotypes
    hammingDistanceBinCounts += HDrowBins

# Make dataframe for distribution of pairwise hamming distances
hammingDistanceBinCounts = np.trim_zeros(hammingDistanceBinCounts, trim='b')
hamDistDF = pd.DataFrame(hammingDistanceBinCounts, columns=['sequence_pairs_with_n_hamming_distance'])
hamDistDF = hamDistDF.reset_index().rename(columns={'index':'n'})
hamDistDF.to_csv(snakemake.output.HamDistCSV, index=False)

# plot hamming distance distribution
plotTitle = f"{tag}_{bc}"
xLabel, yLabel = hamDistDF.columns
TOOLTIPS = [(xLabel, f'@{xLabel}'), (yLabel, f'@{yLabel}')]
hamDistPlot = figure(title=plotTitle, plot_width=600, plot_height=400, x_range=(-0.7,config['hamming_distance_distribution_plot_x_max']), tooltips=TOOLTIPS)
hamDistPlot.vbar(x=xLabel, top=yLabel, width=0.5, bottom=0, color='black', source=ColumnDataSource(hamDistDF))
hamDistPlot.xaxis.axis_label = xLabel
hamDistPlot.yaxis.axis_label = yLabel
output_file(snakemake.output.HamDistPlot)
save(hamDistPlot)

### generate network graph of sequences as nodes and edges connecting all nodes, with inverse hamming distance as edge weight. Plot with holoviews
hammingDistance2Darray = np.vstack(hammingDistanceRows)
hammingDistanceMatrixDF = pd.DataFrame(hammingDistance2Darray, columns=list(genotypesDF.index), index=list(genotypesDF.index))
hammingDistanceMatrixDF.replace(to_replace=-1, value=pd.NA, inplace=True)
hammingDistanceEdgesDF = hammingDistanceMatrixDF.stack().reset_index()
hammingDistanceEdgesDF.columns = ['source', 'target', 'hammingDistance']
hammingDistanceEdgesDF = hammingDistanceEdgesDF[hammingDistanceEdgesDF['hammingDistance'] < max(3, np.argmax(hammingDistanceBinCounts))]    # filter out edges with hamming distance greater than or equal to (a) the maximum hamming distance bincount (will be median for normal distribution), or (b) 3, whichever is larger

def mutCountHDweighting(source,target, hammingDistance):
    sourceMutCount = int(genotypesDF.loc[source,'NT_substitutions_count'])
    targetMutCount = int(genotypesDF.loc[target,'NT_substitutions_count'])
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
nodeColorMap = 'blues' if config['force_directed_plot_node_color']!='barcode(s)' else 'rainbow'

# linear equation to scale Node sizes to range from 5-25 (for node weight)
maxCount = genotypesDF[config['force_directed_plot_node_size']].max()
minCount = genotypesDF[config['force_directed_plot_node_size']].min()
slopeN = (25-5) / (maxCount-minCount)
interceptN = 5 - (slopeN*minCount)

# equation to scale log of edge widths 0.05-7 (for edge weight)
maxWeight = np.log(hammingDistanceEdgesDF['weight'].max())
minWeight = np.log(hammingDistanceEdgesDF['weight'].min())
slopeW = (10-0.05) / (maxWeight-minWeight)
interceptW = 0.1 - (slopeW*minWeight)

# equation to scale log of edge widths 0.05-0.4 (for edge Alpha, or opacity)
slopeA = (0.4-0.05) / np.absolute(maxWeight-minWeight)
interceptA = 0.1 - (slopeA*minWeight)

networkPlot = hv.Graph.from_networkx(G, nx.layout.fruchterman_reingold_layout).opts(
    hv.opts.Graph(node_size=(hv.dim(config['force_directed_plot_node_size'])*slopeN)+interceptN, node_color=config['force_directed_plot_node_color'], cmap=nodeColorMap,
                    edge_line_width=(np.log(hv.dim('weight'))*slopeW)+interceptW, edge_color=np.log(hv.dim('weight')), edge_cmap='Inferno', edge_alpha=np.log(hv.dim('weight'))*slopeA+interceptA))
hv.save(networkPlot, snakemake.output.GraphPlot, backend='bokeh')