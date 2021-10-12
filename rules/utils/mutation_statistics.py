""" script for maple pipeline

Uses output data files from rule mutation_analysis for all files being processed, calculates
interesting statistics from these data files, and outputs these statistics into a .csv file
"""

import numpy as np
import pandas as pd
from Bio import SeqIO
import statistics
import collections

### Asign variables from config file and inputs
config = snakemake.config
inputList = snakemake.input
###


def main():
    statsList = [] # list to be populated with one row of mutation data per tag/barcode combination
    fDict = inFileDict(inputList)
    datatypes = ['genotypes', 'failures', 'NT-muts-frequencies', 'NT-muts-distribution']
    if config['do_AA_analysis']:
        cols = ['tag', 'barcode_group', 'total_seqs', 'total_failed_seqs', 'total_AA_mutations', 'unique_AA_mutations', 'mean_AA_mutations_per_seq', 'median_AA_mutations_per_seq',
        'total_NT_mutations', 'unique_NT_mutations', 'mean_NT_mutations_per_base', 'mean_NT_mutations_per_seq', 'median_NT_mutations_per_seq', 'total_transversions', 'total_transitions', 'unique_transversions', 'unique_transitions']
        datatypes.extend(['AA-muts-distribution', 'AA-muts-frequencies'])
    else:
        cols = ['tag', 'barcode_group', 'total_seqs', 'total_failed_seqs',
        'total_NT_mutations', 'unique_NT_mutations', 'mean_NT_mutations_per_base', 'mean_NT_mutations_per_seq', 'median_NT_mutations_per_seq', 'total_transversions', 'total_transitions', 'unique_transversions', 'unique_transitions']
    for tag in fDict:
        refSeqfasta = config['runs'][tag]['reference']
        referenceLength = len(list(SeqIO.parse(refSeqfasta, 'fasta'))[1].seq)
        for bcGroup in fDict[tag]:
            DFdict = {}
            for dType in datatypes:
                DFdict[dType] = pd.read_csv(fDict[tag][bcGroup][dType], index_col=0)
            
            NTdist = DFdict['NT-muts-distribution']['seqs_with_n_NTsubstitutions']
            totalSeqs = NTdist.sum()
            
            failCount = len(DFdict['failures'])

            
            NTmuts = DFdict['NT-muts-frequencies'].transpose()
            if not config['mutations_frequencies_raw']:
                NTmuts = NTmuts * totalSeqs
                NTmuts = np.rint(NTmuts)
            NTmuts_unique = NTmuts.where(NTmuts == 0, 1) # generate a dataframe that uses only 1s instead of tracking the number of occurences of a particular mutation
            total_NT_mutations = NTmuts.values.sum()
            unique_NT_mutations = NTmuts_unique.values.sum()
            NTmuts.reset_index(inplace=True)
            NTmuts_unique.reset_index(inplace=True)
            transNTmuts = NTmuts.apply(lambda row:
                transversions_transitions(row['index'], row['A'], row['T'], row['G'], row['C']), axis=1, result_type='expand')
            transNTmuts_unique = NTmuts_unique.apply(lambda row:
                transversions_transitions(row['index'], row['A'], row['T'], row['G'], row['C']), axis=1, result_type='expand')
            allMutTypes = NTmuts.apply(lambda row:
                mut_type(row['index'], row['A'], row['T'], row['G'], row['C']), axis=1, result_type='expand')
            allMutTypes_unique = NTmuts_unique.apply(lambda row:
                mut_type(row['index'], row['A'], row['T'], row['G'], row['C']), axis=1, result_type='expand').add_suffix('_unique')
            NTmuts = pd.concat([NTmuts, transNTmuts, allMutTypes], axis='columns')
            NTmuts_unique = pd.concat([NTmuts_unique, transNTmuts_unique, allMutTypes_unique], axis='columns')

            valuesList = [tag, bcGroup, totalSeqs, failCount]

            if config['do_AA_analysis']:
                AAdist = DFdict['AA-muts-distribution']['seqs_with_n_AAsubstitutions']
                if config['mutations_frequencies_raw']:
                    total_AA_mutations = DFdict['AA-muts-frequencies'].values.sum()
                else:
                    total_AA_mutations = round((DFdict['AA-muts-frequencies'] * totalSeqs).values.sum())
                unique_AA_mutations = DFdict['AA-muts-frequencies'].where(DFdict['AA-muts-frequencies'] == 0, 1).values.sum()
                valuesList.extend([total_AA_mutations, unique_AA_mutations, compute_mean_from_dist(AAdist), compute_median_from_dist(AAdist)])

            mean_NT_muts_per_seq = compute_mean_from_dist(NTdist)
            valuesList.extend([total_NT_mutations,  unique_NT_mutations, mean_NT_muts_per_seq/referenceLength, mean_NT_muts_per_seq, compute_median_from_dist(NTdist),
                NTmuts['transversions'].sum(), NTmuts['transitions'].sum(), NTmuts_unique['transversions'].sum(), NTmuts_unique['transitions'].sum()] + [allMutTypes[mutType].sum() for mutType in allMutTypes] + [allMutTypes_unique[mutType].sum() for mutType in allMutTypes_unique])
            
            statsList.append(valuesList)

    cols.extend([column for column in allMutTypes]+[column for column in allMutTypes_unique])
    statsDF = pd.DataFrame(statsList, columns=cols)
    statsDF.sort_values('barcode_group', inplace=True)

    statsDF.to_csv(str(snakemake.output), index=False)

def compute_mean_from_dist(dist):
    """compute mean from pandas series distribution"""
    total = 0
    for n, count in enumerate(dist):
        total += n*count
    if total!=0:
        return round(total/dist.sum(), 2)
    else:
        return 0

def compute_median_from_dist(dist):
    """bad way to compute median from distribution file"""
    seqList = []
    for n, count in enumerate(dist):
        for _ in range(0,count):
            seqList.append(n)
    if all([count==0 for count in dist]):
        return 0
    else:
        return int(statistics.median(seqList))

def mut_type(WT, A, T, G, C):
    """ returns number of each type of mutation as columns
    Used on one sequence position at a time so only one of the four wtNT
    will not be 0 for an individual function call, but combining all outputs
    for all sequence positions gives the total number of each type
    """
    wtNT = WT[0]
    mutsDict = {'A':{}, 'T':{}, 'G':{}, 'C':{}} #nested dict to be used for tracking all 12 types of substitutions
    nts = 'ATGC'
    for wt in nts:
        for mut in nts:
            mutsDict[wt][mut] = 0
    mutsDict[wtNT]['A'] = A
    mutsDict[wtNT]['T'] = T
    mutsDict[wtNT]['G'] = G
    mutsDict[wtNT]['C'] = C
    outDict = collections.OrderedDict()
    for nt in nts:
        for mut in nts:
            if mut!=nt:
                outDict[f'{nt}->{mut}'] = mutsDict[nt][mut]
    return outDict


def transversions_transitions(WT, A, T, G, C):
    wtNT = WT[0]
    transversions = 0
    transitions = 0
    if wtNT in ['A', 'G']:
        transitions += A
        transitions += G
        transversions += T
        transversions += C
    elif wtNT in ['T', 'C']:
        transversions += A
        transversions += G
        transitions += T
        transitions += C
    return {'transversions':transversions, 'transitions':transitions}

def inFileDict(inFileList):
    """ generate a nested dictionary of the input files organized by sample and barcode
        in the format: dict[sample][barcodeGroup][dataType]=fileName """
    outDict = {}
    for f in inFileList:
        sample = f.split('_')[-3].split('/')[-1]
        barcodes = f.split('_')[-2]
        dType = f.split('_')[-1].split('.')[0]
        if sample not in outDict.keys():
            outDict[sample] = {}
        if barcodes not in outDict[sample].keys():
            outDict[sample][barcodes] = {}
        outDict[sample][barcodes][dType] = f
    return outDict

if __name__=='__main__':
    main()