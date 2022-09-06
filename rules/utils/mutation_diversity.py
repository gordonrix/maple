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
from Bio.Seq import translate

### Asign variables from config file
config = snakemake.config
tag = snakemake.wildcards.tag
bc = snakemake.wildcards.barcodes
maxHD = config.get('diversity_plot_hamming_distance_edge_limit', False)
###

genotypesDF = pd.read_csv(snakemake.input[0], na_filter=False)
genotypesDF.drop(genotypesDF[genotypesDF['NT_insertions'] != ''].index, inplace=True)
genotypesDF.drop(genotypesDF[genotypesDF['NT_deletions'] != ''].index, inplace=True)
genotypesDF.drop(genotypesDF[genotypesDF['count'] == 0].index, inplace=True) # removes wildtype row if count is 0

if snakemake.params.downsample:
    if snakemake.params.downsample < len(genotypesDF):
        genotypesDF = genotypesDF.sample(n=snakemake.params.downsample, random_state=1)

genotypesDF.drop(['NT_insertions','NT_deletions'], axis=1, inplace=True)
genotypesDF.set_index('genotype_ID', inplace=True)
refSeqfasta = config['runs'][tag]['reference']
NTref = list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq
NTrefLength = len(NTref)
NTs = 'ATGC'
AAs = 'ACDEFGHIKLMNPQRSTVWY*'
if config['do_AA_mutation_analysis'][tag]:
    AAref = translate(list(SeqIO.parse(refSeqfasta, 'fasta'))[2].seq)
    AArefLength = len(AAref)

# construct a 3 dimensional array of all sequences encoded as 2 dimensional arrays, which are the length of the sequence and all possible mutations, and the third dimension
#   is the genotype ID of the sequence in the order given by the genotypes .csv file, with wildtype (an array of 0s) as the first sequence. From this 3d array,
#   a 2d array will be contructed, consisting of the pairwise hamming distances between each pair of sequences

def sequence_array_from_substitutions(substitutions, NTorAA, array):
    """
    from a str-type list of comma separated substitutions and a numpy zero array,
    will flip any 0s corresponding to a mutation to 1 and return the resulting array

    substitutions   str, comma separated substitutions of the form XNY, where X is wt, Y is the mutation, and N is the index
        N must not exceed the 0-dimension length of `zeroArray`, X and Y may be any character within `NTorAA`
    NTorAA          string of all acceptable characters that may be substitutions
    array       numpy zero array. 1-dimension length must be the same as the length of NTorAA
    """
    if substitutions == '':
        return array
    substitutionsList = substitutions.split(', ')
    for mutation in substitutionsList:
        mut = mutation[-1]
        posi = int(mutation[1:-1])-1
        array[posi,NTorAA.find(mut)] += 1
    return array

def seq_3D_array_from_genotypes_list(genotypesList, refSeq, NTorAA):
    """
    from a list of genotypes represented as a list of comma separated substitutions, generate a 3D array representation of
    this list of genotypes such that an array of zeros of dimensions (len(refSeq), len(NTorAA)) represents the WT sequence
    """
    seqLength = len(refSeq)
    numChars = len(NTorAA)
    numGenotypes = len(genotypesList)
    genotypes3DArray = np.zeros((seqLength, numChars, numGenotypes), dtype=int)

    for i, subsList in enumerate(genotypesList):
        zeroArray = np.zeros((seqLength, numChars), dtype=int)
        genotypes3DArray[:,:,i] = sequence_array_from_substitutions(subsList, NTorAA, zeroArray)
    
    return genotypes3DArray

def HDdist_from_genotypes_list(genotypesList, refSeq, NTorAA):
    """
    from a list of genotypes, the wild type sequence, and a list of letter representations of either nucleotides or amino acids,
    will return a DataFrame of bincounts of pairwise hamming distances among all of the sequences
    """

    muts3Darray = seq_3D_array_from_genotypes_list(genotypesList, refSeq, NTorAA)

    matrixPad = 0 # number of NaNs to pad the beginning of each row so that the hamming distance positions match up with the columns. These are essentially placeholders for hamming distances that were already calculated, and can be found on the other half of the N by N square
    hammingDistanceRows = []
    hammingDistanceBinCounts = np.zeros(len(refSeq), dtype=int)    # Hamming distances for each row converted into bin counts and multiplied by # of times sequence is observed in the dataset

    count = 0
    # iterate through each row and generate counts of hamming distances (bincounts) based on hamming distance between the row genotype and
    #   all other genotypes, including itself, to be used for plotting distribution, and also generate a matrix of hamming distances between
    #   the row genotype and all other genotypes, excluding itself, to be used for plotting a network graph
    for row in range(0, len(genotypesList)):
        matrixPad += 1

        hammingDistance3Darray = np.absolute( muts3Darray[:,:,row:] - muts3Darray[:,:,row].reshape(len(refSeq), len(NTorAA), 1) )       # subtract the 3D array representing remaining genotypes from 3D array representing many copies of the genotype for this row, then take absolute value of all elements in this 3D array
        HDrow = np.sum(hammingDistance3Darray, axis=1)                                                                                  # sum across nucleotides/AAs to get 2D array of total mutations at each position and for each genotype
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
        HDrowBins.resize(len(refSeq))                                                                                       # reshape to be the same length as the sequence, which is the maximum possible hamming distance, allowing for concatenation with bincounts for other genotypes
        hammingDistanceBinCounts += HDrowBins
        count+=1
        print(f"{str(np.round_(count/len(genotypesList)*100))[:-2]}% finished with hamming distance calculation",end='\r')
    print('')

    # Make dataframe for distribution of pairwise hamming distances
    hammingDistanceBinCounts = np.trim_zeros(hammingDistanceBinCounts, trim='b')
    HDdist = pd.DataFrame(hammingDistanceBinCounts, columns=['sequence_pairs_with_n_hamming_distance'])
    HDdist = HDdist.reset_index().rename(columns={'index':'n'})

    return HDdist, hammingDistanceRows

ntHDdist, ntHammingDistanceRows = HDdist_from_genotypes_list(genotypesDF['NT_substitutions'].to_list(), NTref, NTs)
ntHDdist.to_csv(snakemake.output.ntHamDistCSV, index=False)

if config['do_AA_mutation_analysis'][tag]:
    aaHDdist, aaHammingDistanceRows = HDdist_from_genotypes_list(genotypesDF['AA_substitutions_nonsynonymous'].to_list(), AAref, AAs)
    aaHDdist.to_csv(snakemake.output.aaHamDistCSV, index=False)

### generate a dataframe that will get passed to the plotting script
if len(ntHammingDistanceRows)==0:
    exit('[ERROR] mutation_diversity failed because there were no sequence pairs to process after preprocessing.')  
hammingDistance2Darray = np.vstack(ntHammingDistanceRows)
hammingDistanceMatrixDF = pd.DataFrame(hammingDistance2Darray, columns=list(genotypesDF.index), index=list(genotypesDF.index))
hammingDistanceMatrixDF.replace(to_replace=-1, value=pd.NA, inplace=True)
hammingDistanceEdgesDF = hammingDistanceMatrixDF.stack().reset_index()
hammingDistanceEdgesDF.columns = ['source', 'target', 'hammingDistance']
hammingDistanceEdgesDF.to_csv(snakemake.output.edges, index=False)
