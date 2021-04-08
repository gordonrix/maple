"""part of nanopype-MACE pipeline, written by Gordon Rix
UMI_consensus.py
"""

# import gzip
# import re
import pandas as pd
import numpy as np
import pysam
import os
import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
# from demux import BarcodeParser

def main():

    ### Asign variables from config file and input
    config = snakemake.config
    tag = snakemake.wildcards.tag
    BAMin = snakemake.input.grouped
    logIn = snakemake.input.log

    FASTQout = snakemake.output[0]
    minimum = config['UMI_consensus_minimum']
    maximum = config['UMI_consensus_maximum']

    cUMIs = UMI_consensus(config['runs'], tag, BAMin, logIn, FASTQout, minimum, maximum)
    cUMIs.generate_consensus_FASTQ_file()

class UMI_consensus:

    def __init__(self, runsConfig, tag, BAMin, logIn, FASTQout, minimum, maximum):
        """
        arguments:

        runsConfig      - snakemake config dictionary for all runs
        tag             - tag for sequences to be UMI extracted, defined in config file
        BAMin           - BAM input file
        logIn           - log file from UMI group
        BAMout          - BAM file containing a consensus sequence for each UMI group with a final count above the threshold imposed in the config file
        threshold       - minimum number of reads in UMI group to be used for consensus generation. see config['UMI_consensus_threshold']
        """
        self.tag = tag
        self.BAMin = BAMin
        self.logIn = logIn
        self.FASTQout = FASTQout
        self.refSeqfasta = runsConfig[tag]['reference']
        self.reference = list(SeqIO.parse(self.refSeqfasta, 'fasta'))[0]
        self.refStr = str(self.reference.seq)
        self.minimum = minimum
        self.maximum = maximum


    def clean_alignment(self, BAMentry):
        """given a pysam.AlignmentFile BAM entry, uses CIGAR string and reference sequence to align
        the query string to match the reference string positions returns these three strings, the list of quality
        scores, a list of all insertions, and a list of all deletions
        """

        refIndex = BAMentry.reference_start
        queryIndex = 0
        refAln = ''
        queryAln = ''
        queryQualities = []
        insertions = [] # list of tuples where first element is index within reference, second element is sequence inserted

        # for alignments that don't extend to the beginning of the reference sequence, pad with '-'
        refAln += self.refStr[:refIndex]
        padLen = refIndex
        queryAln += '-'*padLen
        queryQualities += [1]*padLen

        for cTuple in BAMentry.cigartuples: # build alignment through combining consecutive segments, which are determined by the cigar tuples, in which the first value determines the operation and the second value determines the length

            if cTuple[0] == 0: #no indel
                refSegment = self.refStr[refIndex:refIndex+cTuple[1]]
                querySegment = BAMentry.query_alignment_sequence[queryIndex:queryIndex+cTuple[1]]
                refAln += refSegment
                queryAln += querySegment
                queryQualities += BAMentry.query_alignment_qualities[queryIndex:queryIndex+cTuple[1]]
                refIndex += cTuple[1]
                queryIndex += cTuple[1]

            elif cTuple[0] == 1: #insertion, not added to sequence to maintain alignment
                insertions.append((refIndex, BAMentry.query_alignment_sequence[queryIndex:queryIndex+cTuple[1]]))
                queryIndex += cTuple[1]

            elif cTuple[0] == 2: #deletion, '-' added to sequence to maintain alignment
                refAln += self.refStr[refIndex:refIndex+cTuple[1]]
                queryAln += '-'*cTuple[1]
                queryQualities += [-1]*cTuple[1]
                refIndex += cTuple[1]

        # for alignments that don't extend to the end of the reference sequence, pad with '-'
        refAln += self.refStr[BAMentry.reference_end:]
        padLen = len(self.refStr)-BAMentry.reference_end
        queryAln += '-'*padLen
        queryQualities += [-1]*padLen

        return  refAln, queryAln, queryQualities, insertions

    # taken from https://stackoverflow.com/questions/51318081/numpy-argmax-when-values-are-equal
    # numpy argmax but will return -1 if max value appears more than once
    @staticmethod
    def argmax_iffOnce(a):
        rows = np.where(a == a.max(axis=1)[:, None])[0]
        rows_multiple_max = rows[:-1][rows[:-1] == rows[1:]]
        argmax = a.argmax(axis=1)
        argmax[rows_multiple_max] = -1
        return argmax

    def generate_consensus_FASTQ_record(self, BAMlist):
        """
        given a list of pysam AlignedSegment objects, returns a string that comprises the consensus sequence for the list, and a matching list of quality scores.
        Consensus for non-insertions is determined primarily by unweighted voting (maximum appearances of nucleotide or gap). If there is a tie
        between two characters, the maximum average quality score for each character will be used instead (gaps are manually assigned a quality score
        of -1, and can thus only be chosen on the basis of the first criterion).
        """
        chars = 'ATGC-'

        # arrays for determining consensus sequence. Unweighted is the sum of occurrences of each possible character (chars) at each position, QscoreWeighted is the average quality score of each possible char at each position
        consensusUnweighted_array = np.zeros((len(self.refStr), len(chars)))
        consensusQscoreWeighted_array = np.zeros((len(self.refStr), len(chars)))

        # nested dict where outer key is the position within the sequence, the inner key is the insertion sequence, and the value is the number of times that insertion position+sequence is observed across all BAM entries for the UMI group
        insertionsDict = {}

        first=True
        for BAMentry in BAMlist:
            refAln, queryAln, Qscores, insertions = self.clean_alignment(BAMentry)
            if first:
                self.out.write(refAln+'\n')
                first=False
            self.out.write(queryAln + '    ' + str(insertions) + '\n')
            for i,nt in enumerate(queryAln):
                consensusUnweighted_array[i][chars.find(nt)] += 1
                consensusQscoreWeighted_array[i][chars.find(nt)] += Qscores[i]
            for ins in insertions:
                index = ins[0]
                seq = ins[1]
                if index not in insertionsDict:
                    insertionsDict[index] = {seq:1}
                elif seq not in insertionsDict[index]:
                    insertionsDict[index][seq] = 1
                else:
                    insertionsDict[index][seq] += 1

        # convert to average, ignore the many instances of 0/0
        with np.errstate(divide='ignore', invalid='ignore'):
            consensusQscoreWeighted_array = consensusQscoreWeighted_array / consensusUnweighted_array

        

        consensusSeq = [] # consensus sequence with deletions as '-' but without insertions
        consensusQuals = []

        for i, (maxU, maxQ) in enumerate(zip(self.argmax_iffOnce(consensusUnweighted_array), np.argmax(consensusQscoreWeighted_array, axis=1))):
            if maxU != -1:  # use maximum number of occurrences is there isn't a tie
                consensusSeq.append(chars[maxU])
                consensusQuals.append(consensusQscoreWeighted_array[i][maxU])
            else:           # use maximum average quality score if there is a tie in number of occurrences
                consensusSeq.append(chars[maxQ])
                consensusQuals.append(maxQ)

        averageQscore = np.mean(consensusQuals) # use average quality score for insertions

        consensusInsertions = {} # dict of insertions that are present in enough sequences to be counted. key is index of insertion and value is insertion sequence
        for insIndex in insertionsDict:
            maxOccurrences = max(insertionsDict[insIndex].values())
            if (maxOccurrences >= 0.5*len(BAMlist)) and (len(BAMlist) > 2): # insertion if found in more than 50% of BAM entries when more than 2 BAMentries used for consensus
                consensusInsertions[insIndex] = max(insertionsDict[insIndex], key=insertionsDict[insIndex].get)

        # add insertions to string in reverse order to maintain index validity
        consensusInsertionsIndexesSorted = list(consensusInsertions.keys())
        consensusInsertionsIndexesSorted.sort(reverse=True)
        for insIndex in consensusInsertionsIndexesSorted:
            insSequence = consensusInsertions[insIndex]
            consensusSeq = consensusSeq[:insIndex] + list(consensusInsertions[insIndex]) + consensusSeq[insIndex:]
            consensusQuals = consensusQuals[:insIndex] + [averageQscore]*len(insSequence) + consensusQuals[insIndex:]

        self.out.write(''.join(consensusSeq))
        self.out.write('\n\n')

        # remove deletions
        consensusQuals = [q for i,q in enumerate(consensusQuals) if consensusSeq[i] != '-']
        consensusSeq = ''.join([c for c in consensusSeq if c != '-'])

        return consensusSeq, consensusQuals

    def generate_consensus_FASTQ_file(self):
        """
        determines which UMI groups to generate consensus sequences for based upon read count threshold then,
        for each UMI group, loops through the input BAM file once and finds all the BAM entries belonging
        to that UMI group, generates a consensus sequence, and writes the consensus sequence as a BAM entry
        with CIGAR strings to the output BAM file
        """

        pysam.index(self.BAMin)
        BAMin = pysam.AlignmentFile(self.BAMin, 'rb')
        logDF = pd.read_csv(self.logIn, sep='\t')
        UMI_groups_above_threshold = logDF[logDF['final_umi_count']>=self.minimum][['final_umi', 'final_umi_count', 'unique_id']].sort_values(['final_umi_count', 'unique_id'], ascending=[False,True]).reset_index(drop=True)
        UMI_groups_above_threshold = UMI_groups_above_threshold.drop_duplicates(subset=['unique_id']).reset_index()
        # add a column to batch sequences into groups to minimize looping through BAM file to find sequences, without putting the whole BAM file in memory
        UMI_groups_above_threshold['batch'] = UMI_groups_above_threshold.apply( lambda row: int(row['index']/10000), axis=1 )

        with open(self.FASTQout, 'a') as output:

            fqBatch = []
            
            for batch_index in range(0, UMI_groups_above_threshold['batch'].max()+1):
                batch_DF = UMI_groups_above_threshold[UMI_groups_above_threshold['batch']==batch_index]
                batch_UMI_IDs = list(batch_DF['unique_id'])
                UMI_BAMbatchDict = {}

                # dictionary that keeps track of the number of sequences being used for the consensus. If this number goes above the set maximum, the associated UMI group ID will be removed from the list of UMI IDs
                UMI_group_counts = {}
                for UMIgroup in batch_UMI_IDs:
                    UMI_group_counts[UMIgroup] = 0

                for BAMentry in BAMin:
                    ID = BAMentry.get_tag('UG')
                    if ID in batch_UMI_IDs:
                        if ID not in UMI_BAMbatchDict:
                            UMI_BAMbatchDict[ID] = []
                        UMI_BAMbatchDict[ID].append(BAMentry)
                        UMI_group_counts[ID] += 1
                        if UMI_group_counts[ID] >= self.maximum:
                            while ID in batch_UMI_IDs: batch_UMI_IDs.remove(ID)
                    if len(batch_UMI_IDs) == 0: # stop searching for more reads if user defined maximum reads have been found for all UMI groups 
                        break
                BAMin.reset()
                
                for ID, UMIgroupBAMentries in UMI_BAMbatchDict.items():
                    self.out = open('consensusDebug.txt', 'a')
                    row = UMI_groups_above_threshold[UMI_groups_above_threshold['unique_id']==ID].iloc[0]
                    self.out.write(f"UMIgroup={ID} sampleid={self.tag} readcount={row['final_umi_count']} umi={row['final_umi']}\n")
                    consensusSeq, consensusQuals = self.generate_consensus_FASTQ_record(UMIgroupBAMentries)
                    record = SeqRecord(Seq(consensusSeq), id=f"UMI_{ID}_consensus_sequence", description=f"sampleid={self.tag} readcount={row['final_umi_count']} umi={row['final_umi']}")
                    record.letter_annotations['phred_quality'] = consensusQuals
                    fqBatch.append(record)

                SeqIO.write(fqBatch, output, 'fastq')
                fqBatch = []

if __name__ == '__main__':
    main()