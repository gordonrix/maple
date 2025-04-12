#!/usr/bin/env python3
"""
script from maple pipeline, written by Gordon Rix
splits a FASTQ file (with unique_id and strand in each header) into
many smaller FASTA files, with records grouped by unique_id.
"""

import pandas as pd
import numpy as np
import os
import shutil
import bisect
from Bio import SeqIO
import sys
import argparse

# Use snakemake variables if available; otherwise, use argparse.
if "snakemake" in globals():
    input_fastq = snakemake.input.fastq
    logIn = snakemake.input.log
    out_dir = os.path.split(snakemake.output.fastas[0])[0]
    minimum = snakemake.params.minimum
    maximum = snakemake.params.maximum
    batches = snakemake.params.batches
    tag = snakemake.wildcards.tag
else:
    parser = argparse.ArgumentParser(description="Split FASTQ by UMI groups into batches for consensus")
    parser.add_argument("--input_fastq", required=True, help="Input FASTQ file (with unique_id in header)")
    parser.add_argument("--log", required=True, help="CSV log file from UMI grouping")
    parser.add_argument("--output_dir", required=True, help="Output directory for batch FASTA files. NOTE: this will be deleted and recreated.")
    parser.add_argument("--minimum", type=int, required=True, help="Minimum number of reads per UMI group")
    parser.add_argument("--maximum", type=int, required=True, help="Maximum number of reads per UMI group")
    parser.add_argument("--batches", type=int, required=True, help="Number of batches to split into")
    args = parser.parse_args()
    input_fastq = args.input_fastq
    logIn = args.log
    out_dir = args.output_dir
    minimum = args.minimum
    maximum = args.maximum
    batches = args.batches

def get_unique_id_and_strand(record):
    unique_id = None
    strand = None
    for token in record.description.split():
        if token.startswith("unique_id="):
            unique_id = token.split("=", 1)[1]
        elif token.startswith("strand="):
            strand = token.split("=", 1)[1]
    if unique_id is None:
        sys.exit("Error: unique_id not found in record " + record.id)
    if strand is None:
        sys.exit("Error: strand not found in record " + record.id)
    return unique_id, strand

def safe_clear_directory(directory, ext=".fasta"):
    """
    Checks if the given directory exists and, if so, verifies that every file (ignoring subdirectories)
    in it ends with the specified extension (default is ".fasta").
    If the check passes, the directory is deleted and re-created.
    If any file does not have the proper extension, the function exits with an error.
    """
    if os.path.exists(directory):
        for item in os.listdir(directory):
            path = os.path.join(directory, item)
            if os.path.isfile(path) and not item.lower().endswith(ext):
                sys.exit(f"Error: {directory} contains files that are not {ext} files. Aborting deletion.")
        shutil.rmtree(directory)
    os.mkdir(directory)

class UMISplitter:
    def __init__(self, FASTQin, logIn, out_dir, minimum, maximum, batches):
        self.FASTQin = FASTQin
        self.logIn = logIn
        self.out_dir = out_dir
        self.minimum = minimum
        self.maximum = maximum
        self.batches = batches

    def split(self):
        logDF = pd.read_csv(self.logIn, sep=",")
        UMI_groups = logDF[logDF['final_umi_count']>=self.minimum][['final_umi_count','unique_id']]
        UMI_groups = UMI_groups.sort_values(['final_umi_count','unique_id'], ascending=[False,True]).reset_index(drop=True)
        if len(UMI_groups) == 0:
            sys.exit("No UMI groups above threshold. Check quality.")
        elif len(UMI_groups) < 1000:
            print('[NOTICE] Fewer than 1000 groups above threshold.')
        UMI_groups = UMI_groups.drop_duplicates(subset=['unique_id']).reset_index()
        # add a column to batch sequences into groups to minimize looping through BAM file to find sequences without putting the whole file in memory
        batchSize = 100000
        UMI_groups['batch'] = UMI_groups.apply(lambda row: int(row['index']/batchSize), axis=1)
        
        safe_clear_directory(self.out_dir, ext=".fasta")
        split_fasta_dict = {}
        for x in range(0, self.batches):
            fastaOutName = f'{self.out_dir}/batch{x}.fasta'
            split_fasta_dict[x] = open(fastaOutName, 'w')

        # loop through the BAM file once per batch of UMI IDs
        for batch_index in range(0, UMI_groups['batch'].max()+1):
            batch_DF = UMI_groups[UMI_groups['batch']==batch_index]
            batch_UMI_IDs = set(batch_DF['unique_id'].astype(str))
            # For each batch, build a dictionary keyed by unique_id to hold records and their quality.
            group_rec_dict = {}     
            group_strand_dict = {}      # dict to keep track of how many fwd/rvs strands have been encountered to aim for similar amounts to reduce systematic errors from strand bias
            group_qual_dict = {}        # dict to keep track of average quality scores of reads

            with open(self.FASTQin, "r") as fh:
                for record in SeqIO.parse(fh, "fastq"):
                    uid, strand = get_unique_id_and_strand(record)
                    if uid in batch_UMI_IDs:
                        strand = 1 if strand == '+' else -1
                        mean_q = np.mean(record.letter_annotations["phred_quality"])
                        if uid not in group_rec_dict:
                            group_rec_dict[uid] = []
                            group_strand_dict[uid] = []
                            group_qual_dict[uid] = []
                            
                        # Insert record in sorted order (ascending by quality)
                        idx = bisect.bisect_left(group_qual_dict[uid], mean_q)
                        group_rec_dict[uid].insert(idx, record)
                        group_strand_dict[uid].insert(idx, strand)
                        group_qual_dict[uid].insert(idx, mean_q)
                        
                        # If the group exceeds the maximum, remove reads to increase average quality and (if possible) decrease strand bias
                        if len(group_rec_dict[uid]) > self.maximum:
                            strand_bias = sum(group_strand_dict[uid])
                            if len(set(group_strand_dict[uid])) == 2:
                                if strand_bias < 0:
                                    remove_idx = group_strand_dict[uid].index(-1)
                                elif strand_bias >= 0:
                                    remove_idx = group_strand_dict[uid].index(1)
                            else:
                                remove_idx = 0
                            group_rec_dict[uid].pop(remove_idx)
                            group_strand_dict[uid].pop(remove_idx)
                            group_qual_dict[uid].pop(remove_idx)

            # Write out the records for each unique_id in this batch.
            for uid, rec_list in group_rec_dict.items():
                out_batch = int(uid) % self.batches
                for rec in rec_list:
                    split_fasta_dict[out_batch].write(f'>UMI-{uid}_{rec.id}\n{rec.seq}\n')
        
        for handle in split_fasta_dict.values():
            handle.close()

def main():
    splitter = UMISplitter(input_fastq, logIn, out_dir, minimum, maximum, batches)
    splitter.split()

if __name__ == '__main__':
    main()
