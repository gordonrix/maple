"""
script from maple pipeline, written by Gordon Rix
splits BAM file into many smaller files, which are each separately used to generate consensus sequence,
which are then combined to generate a single consensus.fasta file

"""

# import gzip
# import re
import pandas as pd
import numpy as np
import pysam
import os
import shutil
import datetime
import bisect
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def main():

    ### Asign variables from config file and input
    tag = snakemake.wildcards.tag
    BAMin = snakemake.input.grouped
    logIn = snakemake.input.log

    outDir = os.path.split(snakemake.output.fastas[0])[0]
    minimum = snakemake.params.minimum
    maximum = snakemake.params.maximum
    batches = snakemake.params.batches

    BAMs = UMIBAMs(tag, BAMin, logIn, outDir, minimum, maximum, batches)
    BAMs.split()

class UMIBAMs:

    def __init__(self, tag, BAMin, logIn, outDir, minimum, maximum, batches):
        """
        arguments:

        tag             - tag for sequences to be UMI extracted, defined in config file
        BAMin           - BAM input file
        logIn           - log file from UMI group
        FASTAoutDir     - Directory in which to place all of the fasta file outputs. Each output fasta corresponds to one group of fasta sequences,
                            which are grouped according to the UMI identifiedin the sequence
        minimum         - minimum number of reads a UMI group must have to be used for consensus generation.
        maximum         - maximum number of reads in UMI group to be used for consensus generation. see config['UMI_consensus_maximum']
        batches         - number of files to split up BAM file into, based on modulo of UMI ID (BAM tag 'UG')
        """
        self.tag = tag
        self.BAMin = BAMin
        self.logIn = logIn
        self.tempDir = outDir
        self.minimum = minimum
        self.maximum = maximum
        self.batches = batches

    def split(self):
        """
        determines which UMI groups to generate consensus sequences for based upon minimum and maximum set read count then,
        for each UMI group, loops through the input BAM file and finds all the BAM entries belonging
        to that UMI group, and writes the reads as a BAM file, then sorts the BAM file according to the 'UG' tag
        """

        logDF = pd.read_csv(self.logIn, sep='\t')
        UMI_groups_above_threshold = logDF[logDF['final_umi_count']>=self.minimum][['final_umi', 'final_umi_count', 'unique_id']].sort_values(['final_umi_count', 'unique_id'], ascending=[False,True]).reset_index(drop=True)

        if len(UMI_groups_above_threshold) == 0:
            raise RuntimeError(f"No reads with UMI counts above UMI threshold. Threshold may be too high or sequencing run was of poor quality. Examine `plots/{self.tag}_UMIgroup-distribution` and plots in `plots/nanoplot/` directory to determine the root of the problem.")
        elif len(UMI_groups_above_threshold) < 1000:
            print('[WARNING] Fewer than 1000 reads with UMI counts above UMI threshold. Threshold may be too high or sequencing run was of poor quality. Examine `plots/{self.tag}_UMIgroup-distribution` and plots in `plots/nanoplot/` directory to determine if there is a problem.')

        UMI_groups_above_threshold = UMI_groups_above_threshold.drop_duplicates(subset=['unique_id']).reset_index()
        # add a column to batch sequences into groups to minimize looping through BAM file to find sequences, without putting the whole BAM file in memory
        batchSize = 100000
        UMI_groups_above_threshold['batch'] = UMI_groups_above_threshold.apply( lambda row: int(row['index']/batchSize), axis=1 )
        
        shutil.rmtree(self.tempDir, ignore_errors=True)
        os.mkdir(self.tempDir)

        with pysam.AlignmentFile(self.BAMin, 'rb') as BAMin:

            # open all BAM output files for writing
            splitFastaDict = {}
            for x in range(0, self.batches):
                fastaOutName = f'{self.tempDir}/batch{x}.fasta'
                splitFastaDict[x] = open(fastaOutName, 'w')

            # loop through the BAM file once per batch of sequences. note that this is not the same as the batches of output files
            for batch_index in range(0, UMI_groups_above_threshold['batch'].max()+1):

                batch_DF = UMI_groups_above_threshold[UMI_groups_above_threshold['batch']==batch_index]
                batch_UMI_IDs = list(batch_DF['unique_id'])
                UMI_BAMbatchDict = {}
                UMI_strandTrackDict = {}     # dict to keep track of how many fwd/rvs strands have been encountered to aim for similar amounts to reduce systematic errors from strand bias. fwd recorded as 1, rvs recorded as -1 such that the sum of the list indicates the bias
                UMI_qualityTrackDict = {}    # dict to keep track of average quality scores of reads

                for BAMentry in BAMin:
                    ID = BAMentry.get_tag('UG')
                    if ID in batch_UMI_IDs:
                        if ID not in UMI_BAMbatchDict:
                            UMI_BAMbatchDict[ID] = []
                            UMI_strandTrackDict[ID] = []
                            UMI_qualityTrackDict[ID] = []
                        if BAMentry.query_qualities is not None:
                            BAMmeanQscore = np.mean(BAMentry.query_qualities)
                            insertIndex = bisect.bisect_left(UMI_qualityTrackDict[ID], BAMmeanQscore) # grow list of BAM entries, Qscores, and strand tracking in order of Qscores so that lowest Qscore sequences get removed first
                        else:
                            insertIndex, BAMmeanQscore = 0, 0
                        UMI_BAMbatchDict[ID].insert(insertIndex, BAMentry)
                        UMI_qualityTrackDict[ID].insert(insertIndex, BAMmeanQscore)
                        if BAMentry.is_reverse:
                            UMI_strandTrackDict[ID].insert(insertIndex, -1)
                        else:
                            UMI_strandTrackDict[ID].insert(insertIndex, 1)

                        # UMI IDs will be removed from list if the maximum # of reads have been reached and the difference
                            # in counts of reads of one strand vs another is either 0 (for even # of maximum reads) or 1 (for odd # of maximum reads)
                            # As more reads over the maximum are encountered, reads from the left are removed to approach a strand bias of 0 and to increase the minimum average read quality score
                        UMI_group_strandBias = sum(UMI_strandTrackDict[ID])
                        if len(UMI_BAMbatchDict[ID]) > self.maximum:
                            if UMI_group_strandBias < 0:
                                removeIndex = UMI_strandTrackDict[ID].index(-1) # index searches from the left so lowest quality score with desired strandedness will be removed first
                                UMI_BAMbatchDict[ID].pop(removeIndex)
                                UMI_strandTrackDict[ID].pop(removeIndex)
                                UMI_qualityTrackDict[ID].pop(removeIndex)
                            elif UMI_group_strandBias >= 0:
                                removeIndex = UMI_strandTrackDict[ID].index(1)
                                UMI_BAMbatchDict[ID].pop(removeIndex)
                                UMI_strandTrackDict[ID].pop(removeIndex)
                                UMI_qualityTrackDict[ID].pop(removeIndex)
                BAMin.reset()  # allows for looping through again for the next batch

                for ID, UMIgroupBAMentries in UMI_BAMbatchDict.items():
                    x = int(ID)%self.batches
                    for BAMentry in UMIgroupBAMentries:
                        splitFastaDict[x].write(f'>UMI-{ID}_{BAMentry.qname}\n{BAMentry.query_sequence}\n')
                
        # close output files
        for key in splitFastaDict:
            splitFastaDict[key].close()

    def run_maple_smolecule(self, consensusOutFasta, model, threads, flags):
        """
        loops through all of the BAM files that have been generated by split method and runs a modified medaka smolecule script that accepts
        a BAM file instead of a fasta files, then combines the .fasta outputs and compresses them

        arguments:

        consensusOutFasta   - str, fasta file name to combine all output consensus fasta files into
        model               - medaka model to use for consensus calling
        threads             - number of threads to use
        flags               - flags to use for medaka smolecule
        """
        with open(consensusOutFasta, 'w') as consensuses:
            for key in self.splitBAMdict.keys():
                BAMfile = f'{self.tempDir}/batch{key}_sorted.bam'
                os.system(f"medaka maple_smolecule {self.tempDir} {BAMfile} --model {model} --threads {threads} {flags}")
                batchConsensuses = os.path.join(self.tempDir, 'consensus.fasta')
                with open(batchConsensuses, 'r') as batchConsHandle:
                    consensuses.write(''.join(batchConsHandle.readlines()))
                os.remove(os.path.join(self.tempDir, 'consensus.hdf')) # stale consensus.hdf can sometimes cause problems
                
        os.system(f'gzip {consensusOutFasta}')

if __name__ == '__main__':
    main()