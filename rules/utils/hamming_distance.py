#
#  DESCRIPTION   : Custom python script used by the Maple snakemake pipeline.
#                   Uses vectorized hamming distance calculation and bincounts
#                   to generate a matrix of all pairwise hamming distances and
#                   the distribution of all pairwise hamming distance counts
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

import pandas as pd
import numpy as np
import sklearn
from Bio import SeqIO
from Bio.Seq import translate
from dimension_reduction_genotypes import SequenceEncoder
from dimension_reduction_genotypes import seq_array_from_genotypes

def main():

    genotypesDF = pd.read_csv(snakemake.input[0])

    if snakemake.wildcards.NTorAA == 'NT':
        refSeq = list(SeqIO.parse(snakemake.params.refSeqs, 'fasta'))[1].seq
    elif snakemake.wildcards.NTorAA =='AA':
        refSeq = translate(list(SeqIO.parse(snakemake.params.refSeqs, 'fasta'))[2].seq)
    matrixDF, HDdistDF = HD_matrix_and_dist(refSeq, genotypesDF, snakemake.wildcards.NTorAA, snakemake.params.downsample)
    matrixDF.to_csv(snakemake.output.HDmatrixCSV, index=False)
    HDdistDF.to_csv(snakemake.output.HDdistCSV, index=False)

def pairwise_hamming_distance_matrix(seqArray):
    """
    Given an array of shape (N,L) containing N integer-encoded sequences of length L,
    computes the pairwise hamming distance of all pairs of sequences and outputs these as
    a matrix of shape (N,N) containing all pairwise hamming distances
    """
    # sklearn pairwise distances returns distances as a fraction of the maximum distance,
    #   so multiplying by the maximum distance, then rounding, and converting to int
    return np.rint(sklearn.metrics.pairwise_distances(seqArray,metric='hamming')*seqArray.shape[1]).astype(int)


def bincount_2D(matrix):
    """
    matrix:     2d array

    given a matrix (2d array) of positive integer values of shape (M,N),
        returns bincounts for the matrix as a 2d array of shape
        (M, matrix.max()+1). Essentially a vectorized version of
        np.bincount applied to each element along the 0th dimension of a 2d array

    taken from https://stackoverflow.com/questions/40591754/vectorizing-numpy-bincount
    """
    M = matrix.shape[0]
    maxVal = matrix.max()+1
    matrix = matrix.transpose()

    arr = matrix + (maxVal*np.arange(M))
    return np.bincount(arr.ravel(), minlength=M*maxVal).reshape(M,-1)


def HD_matrix_and_dist(refSeq, genotypesDF, NTorAA, downsample):
    """
    refSeq:         string, reference sequence used to analyze mutations
    genotypesDF:    pandas.DataFrame, genotypes table from the genotypes.csv file
    NTorAA:         string, 'NT' or 'AA', determines whether to calculate
                        NT hamming distance or AA hamming distance
    downsample:     False or int, if int, after removal of genotypes with insertions or deletions,
                        the number of sequences will be reduced down to this number
    """

    seqArray, genotypesUsed = seq_array_from_genotypes(refSeq, genotypesDF, NTorAA)

    if downsample:
        if downsample < len(seqArray):
            idx = np.sort(np.random.choice(np.arange(len(seqArray)), size=downsample, replace=False))
            idx=np.arange(6)
            seqArray = seqArray[idx,:]
            genotypesUsed = genotypesUsed.iloc[idx]

    HDmatrix = pairwise_hamming_distance_matrix(seqArray)
    matrixDF = pd.DataFrame(HDmatrix, index=genotypesUsed.index, columns=genotypesUsed.index).reset_index().rename(columns={'index':'genotype_ID'})

    triangle = np.tril(HDmatrix)            # matrix is symmetric along diagonal, so zero out hamming distances from upper triangle
    HDbincount = bincount_2D(triangle)
    HDbincount[:,0] = np.subtract(HDbincount[:,0], np.arange(HDbincount.shape[0],0,-1))     # remove 0 counts resulting from zeroing out the upper triangle
    counts = genotypesUsed['count'].to_numpy().reshape(len(genotypesUsed),1)
    HDbincount = np.multiply(HDbincount, counts)                                            # multiply hamming distance bincounts for each sequence by the counts for each sequence
    maxHD = HDbincount.shape[1]
    HDs = np.arange(maxHD).reshape(1,maxHD)
    HDdist = HDbincount.sum(0).reshape(1,maxHD)
    HDproportion = np.divide(HDdist, HDdist.sum())
    HDdistDF = pd.DataFrame( np.concatenate((HDs, HDdist, HDproportion), axis=0).T, columns = ['n', 'number_of_sequence_pairs_with_n_hamming_distance', 'proportion_of_sequence_pairs_with_n_hamming_distance'])
    
    return matrixDF, HDdistDF

if __name__ == '__main__':
    main()