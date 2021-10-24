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
genotypesDF.drop(['NT_insertions','NT_deletions'], axis=1, inplace=True)
genotypesDF.set_index('genotype', inplace=True)
refSeqfasta = config['runs'][tag]['reference']
NTref = list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq
NTrefLength = len(NTref)
NTs = 'ATGC'

# construct a 3 dimensional array of all sequences encoded as 2 dimensional arrays, which are the length of the sequence and all four possible nucleotides, and the third dimension
#   is the genotype ID of the sequence in the order given by the genotypes .csv file, with wildtype (an array of 0s) as the first sequence. From this 3d array,
#   a 2d array will be contructed, consisting of the pairwise hamming distances between each pair of sequences
genotypes3DArray = np.zeros((int(len(NTref)), len(NTs), len(genotypesDF)), dtype=int)

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

# posiCount, NTcount, genotypesCount = 0,0,0
# for position in genotypes3DArray[:,0,0]:
#     posiCount += 1
# for NT in genotypes3DArray[0,:,0]:
#     NTcount += 1
# for genotype in genotypes3DArray[0,0,:]:
#     genotypesCount += 1
# # print(posiCount, NTcount, genotypesCount)

hammingDistanceDF = np.empty((len(genotypesDF), len(genotypesDF)), dtype=int) * np.nan

# get hamming distances by finding the difference between arrays for each sequence
for row in range(0, len(genotypesDF)):
    for col in range(0,len(genotypesDF)):
        if col > row: # only calculate once per pair
            hammingDistanceDF[row,col] = np.sum( np.absolute( genotypes3DArray[:,:,row]-genotypes3DArray[:,:,col] ) )

# print(hammingDistanceDF)

matrixPad = 0 # number of NaNs to pad the beginning of each row so that the hamming distance positions match up with the columns. These are essentially placeholders for hamming distances that were already calculated, and can be found on the other half of the N by N square
hammingDistanceRows = []
hammingDistanceBinCounts = np.zeros(NTrefLength, dtype=int)    # Hamming distances for each row converted into bin counts and multiplied by # of times sequence is observed in the dataset

# iterate through each row and generate counts of hamming distances (bincounts) based on hamming distance between the row genotype and
#   all other genotypes, including itself, to be used for plotting distribution, and also generate a matrix of hamming distances between
#   the row genotype and all other genotypes, excluding itself, to be used for plotting a network graph
for row in range(0, len(genotypesDF)):
    matrixPad += 1

    rowGenotype3Darray = np.concatenate([genotypes3DArray[:,:,row][..., np.newaxis]]*(len(genotypesDF)-row), axis=2)    # expand the 2D array of the row genotype to a 3D array that matches the shape of the remaining genotypes
    hammingDistance3Darray = np.absolute( rowGenotype3Darray - genotypes3DArray[:,:,row:] )                             # subtract the 3D array representing remaining genotypes from 3D array representing many copies of the genotype for this row, then take absolute value of all elements in this 3D array
    HDrow = np.sum(hammingDistance3Darray, axis=1)                                                                      # sum across nucleotides to get 2D array of total mutations at each position and for each genotype
    HDrow[HDrow>1] = 1                                                                                                  # set values above 1 to 1. otherwise, genotypes that both have different mutations at a single position will result in a hamming distance of 2 for that position, but hamming distance can only ever be 1 for a single position
    HDrow = np.sum(HDrow, axis=0)                                                                                       # sum across position axis to give 1D array of total hamming distance for each genotype
    hammingDistanceRows.append( np.pad(HDrow[1:], ((matrixPad,0)), constant_values=-1) )                                    # make new row of all hamming distances for the row genotype, padding the beginning of the row with the null value -1 so that when combined the column positions will line up

    genotypeCount = int(genotypesDF.iloc[[row]]['count'])
    HDcountsList = [0] * int((genotypeCount*(genotypeCount-1))/2)                                                          # 0 count for comparison of the row genotype to itself calculated as n(n-1)/2, see https://stackoverflow.com/questions/18859430/how-do-i-get-the-total-number-of-unique-pairs-of-a-set-in-the-database
    # print(genotypesDF.iloc[[row]].index)
    # print('genotypeCount, len(HDcountsList)', genotypeCount, len(HDcountsList))
    # loop through hamming distances for row except first one which is always 0
    for HDindex in range(1, len(HDrow)):

        # increase counts of each hamming distance by the product of the counts for the two sequences being compared
        HDcountsList.extend( [HDrow[HDindex]] * int(genotypesDF.iloc[[row]]['count']) * int(genotypesDF.iloc[[HDindex+row]]['count']) )

    HDrowBins = np.bincount(HDcountsList)                                                                               # get counts of each hamming distance observed, then multiply these bincounts by the count for the row genotype to give the total count of these observed hamming distances
    # print(HDrow)
    # print(HDcountsList)
    # print(HDrowBins)
    # print('')
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
hammingDistanceEdgesDF = hammingDistanceEdgesDF[hammingDistanceEdgesDF['hammingDistance']<=3]    # filter out edges with hamming distance greater than 3

def mutCountHDweighting(source,target, hammingDistance):
    sourceMutCount = len(genotypesDF.loc[source,'NT_substitutions'].split(', '))
    targetMutCount = len(genotypesDF.loc[target,'NT_substitutions'].split(', '))
    return 2**((sourceMutCount+targetMutCount)/2-hammingDistance)

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

maxCount = genotypesDF['count'].max()           # maximum genotype used to scale node sizes. Largest node size is 23
networkPlot = hv.Graph.from_networkx(G, nx.layout.fruchterman_reingold_layout).opts(
    hv.opts.Graph(node_size=3+(hv.dim('count')*(20/maxCount)), edge_color='weight', cmap='Set1',
                    edge_cmap='viridis', edge_line_width=hv.dim('weight')*(15/hammingDistanceEdgesDF['weight'].max())))
hv.save(networkPlot, snakemake.output.GraphPlot, backend='bokeh')