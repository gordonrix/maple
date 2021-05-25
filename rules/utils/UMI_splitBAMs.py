"""part of nanopype-MACE pipeline, written by Gordon Rix
UMI_BAMtoFASTA.py
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

def main():

    ### Asign variables from config file and input
    config = snakemake.config
    tag = snakemake.wildcards.tag
    BAMin = snakemake.input.grouped
    logIn = snakemake.input.log

    outDir = snakemake.output[0]
    minimum = config['UMI_consensus_minimum']
    maximum = config['UMI_consensus_maximum']

    BAMs = splitBAM(tag, BAMin, logIn, outDir, minimum, maximum)
    BAMs.split()

class splitBAM:

    def __init__(self, tag, BAMin, logIn, outDir, minimum, maximum):
        """
        arguments:

        tag             - tag for sequences to be UMI extracted, defined in config file
        BAMin           - BAM input file
        logIn           - log file from UMI group
        FASTAoutDir     - Directory in which to place all of the fasta file outputs. Each output fasta corresponds to all reads for one UMI group
        maximum         - maximum number of reads in UMI group to be used for consensus generation. see config['UMI_consensus_maximum']
        """
        self.tag = tag
        self.BAMin = BAMin
        self.logIn = logIn
        self.outDir = outDir
        self.minimum = minimum
        self.maximum = maximum

    def split(self):
        """
        determines which UMI groups to generate consensus sequences for based upon minimum and maximum set read count then,
        for each UMI group, loops through the input BAM file and finds all the BAM entries belonging
        to that UMI group, and writes the reads as a fasta file
        """

        pysam.index(self.BAMin)
        BAMin = pysam.AlignmentFile(self.BAMin, 'rb')
        logDF = pd.read_csv(self.logIn, sep='\t')
        UMI_groups_above_threshold = logDF[logDF['final_umi_count']>=self.minimum][['final_umi', 'final_umi_count', 'unique_id']].sort_values(['final_umi_count', 'unique_id'], ascending=[False,True]).reset_index(drop=True)

        if len(UMI_groups_above_threshold) == 0:
            raise RuntimeError(f"No reads with UMI counts above UMI threshold. Threshold may be too low or sequencing run was of poor quality. Examine `plots/{self.tag}_UMIgroup-distribution` and plots in `plots/nanoplot/` directory to determine the root of the problem.")
        elif len(UMI_groups_above_threshold) < 1000:
            print('[WARNING] Fewer than 1000 reads with UMI counts above UMI threshold. Threshold may be too low or sequencing run was of poor quality. Examine `plots/{self.tag}_UMIgroup-distribution` and plots in `plots/nanoplot/` directory to determine if there is a problem.')

        UMI_groups_above_threshold = UMI_groups_above_threshold.drop_duplicates(subset=['unique_id']).reset_index()
        # add a column to batch sequences into groups to minimize looping through BAM file to find sequences, without putting the whole BAM file in memory
        batchSize = 100000
        UMI_groups_above_threshold['batch'] = UMI_groups_above_threshold.apply( lambda row: int(row['index']/batchSize), axis=1 )
        
        os.mkdir(self.outDir)
        # loop through the BAM file once per batch of sequences
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
            BAMin.reset()  # allows for looping through again for the next batch
            
            for ID, UMIgroupBAMentries in UMI_BAMbatchDict.items():
                row = UMI_groups_above_threshold[UMI_groups_above_threshold['unique_id']==ID].iloc[0]
                BAMout = pysam.AlignmentFile(os.path.join(self.outDir, f"UMI_{ID}_{row['final_umi']}_reads.bam"), 'wb', template=BAMin)
                for BAMentry in UMIgroupBAMentries:
                    BAMout.write(BAMentry)

if __name__ == '__main__':
    main()