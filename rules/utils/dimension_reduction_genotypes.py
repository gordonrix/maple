#
#  DESCRIPTION   : Custom python script used by the Maple snakemake pipeline.
#                   Encodes sequences as integer arrays, then uses the
#                   PaCMAP dimensionality reduction tool to reduce these 
#                   arrays down to 2 dimensions
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

from logging import exception
import pacmap as pm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from copy import deepcopy
from Bio import SeqIO
import re

def main(ntref, inputCSV, outputCSV):
    """
    ntref       string, nucleotide reference sequence
    inputCSV    string, .csv file name that contains a NT_substitutions, NT_insertions, and NT_deletions columns
    outputCSV   string, .csv file name to output to

    Uses NT substitutions and NT deletions to encode all sequences without insertions as an array, then uses
        PaCMAP to reduce dimensions to 2. These 2 dimensional coordinates are then appended to the original
        .csv file and output as a new .csv file
    """

    genotypes = pd.read_csv(inputCSV)
    array_of_seqs, genotypesNoIns = seq_array_from_genotypes(ntref, genotypes, 'NT')

    # initializing the pacmap instance
    embedding = pm.PaCMAP(n_components=2, MN_ratio=0.5, FP_ratio=2.0) 

    # fit the data to 2 dimensions, add as columns to the original genotypes DF, and export
    seq2D = embedding.fit_transform(array_of_seqs, init="pca")
    seq2Ddf = pd.DataFrame(seq2D, columns=['dim1','dim2'], index=genotypesNoIns.index)
    genotypes = pd.merge(genotypes, seq2Ddf, how='left', left_index=True, right_index=True)
    genotypes.to_csv(outputCSV, index=False)

def seq_array_from_genotypes(refSeq, genotypes, NTorAA):
    """
    args:
        genotypes:          pandas DataFrame of a genotypes.csv file
        refSeq:             string, AA or NT sequence that serves as a template for the mutations
                                in substitutionsList and deletionsList
        NTorAA:             string, either 'NT' or 'AA
    
    returns an array of integer-encoded sequences of shape (N,L) where
        N = len(substitutionsList = len(deletionsList) and L = len(refSeq),
        and also the genotypes DataFrame because it gets modified
    """

    genotypes = genotypes[genotypes['NT_insertions'].isna()]   # remove genotypes with insertions because they mess with indexing and are hard
    
    if NTorAA=='NT':
        chars = 'ATGC-'
        subs = 'NT_substitutions'
    if NTorAA=='AA':
        chars = 'ACDEFGHIKLMNPQRSTVWY*' # no deletions for AA because they are hard
        subs = 'AA_substitutions_nonsynonymous'
        genotypes = genotypes[genotypes['NT_deletions'].isna()]   # remove genotypes with deletions for AA analysis because they are hard

    baseSeq = SequenceEncoder(refSeq,chars)                       # base sequence from reference sequence that will be copied and modified

    # make an array of integer encoded genotypes of shape (N,L),
    #   where N = number of genotypes and L = nucleotide length
    intSeqs = []
    for subs,dels in genotypes[[subs,'NT_deletions']].itertuples(index=False, name=None):
        seq = deepcopy(baseSeq)
        seq.genotype_modify(subs,dels)
        intSeqs.append(seq.integer)

    return np.array(intSeqs), genotypes

class SequenceEncoder:
    def __init__(self, sequence, characters):
        
        #get sequence into an array
        self.seqArray = np.array(list(sequence.upper()))

        # create an encoder dict for consistent encoding
        self.encoderDict = {char:index for index,char in enumerate(characters.upper())}
        
        # create both integer- and onehot-encoded arrays for the sequence and add to self
        self.encode_integer()
        self.encode_onehot()

    def encode_integer(self):
        # integer encode the sequence (eg ATTGC becomes 12234 with characters='ATGC')
        integerEncodedSeq = np.ndarray(self.seqArray.shape, dtype=np.int8)
        for char in self.encoderDict:
            integerEncodedSeq[self.seqArray==char] = self.encoderDict[char]

        self.integer = integerEncodedSeq

    def encode_onehot(self):
        # use integer encoded sequence to generate one hot encoded sequence
        self.onehot = np.zeros((len(self.integer), len(self.encoderDict)), dtype=np.int8)
        self.onehot[np.arange(len(self.integer)), self.integer] = 1
    
    def genotype_modify(self, substitutions, deletions):
        """
        modify the sequence array and integer-encoded sequence according to provided substitutions and deletions,
            then update one hot -encoded sequence
                args:
                    substitutions    string, comma separated list of 1-index mutations, ex: "A234C, T301G"
                                        first and last characters must be present in `characters`
                    deletions        string, comma separated list of 1-index nt deletions, ex: "200del3, 206del1"
        """

        if type(substitutions) == str:
            for mut in substitutions.replace(' ','').split(','):
                posi = int(mut[1:-1])-1 # convert to 0-index
                mutNT = mut[-1]
                self.seqArray[posi] = mutNT
                self.integer[posi] = self.encoderDict[mutNT]
        if type(deletions) == str:
            for deletion in deletions.replace(' ','').split(','):
                posi, number = [int(x) for x in deletion.split('del')]
                posi = posi-1 # convert to 0-index
                self.seqArray[posi:posi+number] = '-'
                self.integer[posi:posi+number] = self.encoderDict['-']
        
        self.encode_onehot()

if __name__ == '__main__':

    ntref = str(list(SeqIO.parse(snakemake.params.refSeqs, 'fasta'))[1].seq)
    main(ntref, snakemake.input.genotypes, snakemake.output.reduced)