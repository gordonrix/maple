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
from Bio import SeqIO
from Bio.Seq import translate
from dimension_reduction_genotypes import sequenceEncoder
from dimension_reduction_genotypes import seq_array_from_genotypes

def pairwise_hamming_distance_matrix(seqArray):
    """
    Given an array of shape (N,L) containing N integer-encoded sequences of length L,
    computes the pairwise hamming distance of all pairs of sequences and outputs these as
    a matrix of shape (N,N) containing all pairwise hamming distances

    taken from https://stackoverflow.com/questions/42752610/python-how-to-generate-the-pairwise-hamming-distance-matrix
    """
    return (seqArray[:, None, :] != seqArray).sum(2)


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
                        NT hamming distance or AA haming distance
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
    HDdist = HDbincount.sum(0).reshape(1,maxHD)
    HDs = np.arange(maxHD).reshape(1,maxHD)
    HDdist = np.concatenate((HDs,HDdist), axis=0).T
    HDdistDF = pd.DataFrame(HDdist, columns = ['n', 'number_of_sequence_pairs_with_n_hamming_distance'])
    
    return matrixDF, HDdistDF


def main():

    genotypesDF = pd.read_csv(snakemake.input[0])

    NTref = list(SeqIO.parse(snakemake.params.refSeqs, 'fasta'))[1].seq
    NT_matrixDF, NT_HDdistDF = HD_matrix_and_dist(NTref, genotypesDF, 'NT', snakemake.params.downsample)
    NT_matrixDF.to_csv(snakemake.output.ntHDmatrixCSV, index=False)
    NT_HDdistDF.to_csv(snakemake.output.ntHDdistCSV, index=False)

    if snakemake.params.do_AA:
        AAref = translate(list(SeqIO.parse(snakemake.params.refSeqs, 'fasta'))[2].seq)
        AA_matrixDF, AA_HDdistDF = HD_matrix_and_dist(AAref, genotypesDF, 'AA', snakemake.params.downsample)
        AA_matrixDF.to_csv(snakemake.output.aaHDmatrixCSV, index=False)
        AA_HDdistDF.to_csv(snakemake.output.aaHDdistCSV, index=False)

if __name__ == '__main__':
    main()