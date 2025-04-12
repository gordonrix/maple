#!/usr/bin/env python3
"""
transfer_cluster_id.py

Reads two FASTQ files simultaneously:
  - One FASTQ file with UMI information.
  - One FASTQ file with full-length sequences.
For each pair of records, it extracts the 'cluster_id' from the UMI record's description
(which must contain a field like "cluster_id=XXXX") and appends that cluster_id to the description
of the corresponding sequence record. cluster_id is added as a new field "unique_id=XXXX".
If any UMI record is missing a cluster_id or if the two files have different numbers of records,
the script exits with an error.
The modified sequence records are then written to an output FASTQ file.
"""

import sys
from itertools import zip_longest
from Bio import SeqIO
import argparse

# If running under Snakemake, the global variable 'snakemake' will be available.
if 'snakemake' in globals():
    class Args: pass
    args = Args()
    args.umi_fastq = snakemake.input.umi_fastq
    args.seq_fastq = snakemake.input.seq_fastq
    args.output_fastq = snakemake.output.fastq
else:
    parser = argparse.ArgumentParser(description="Transfer cluster_id from UMI FASTQ to sequence FASTQ")
    parser.add_argument("--umi_fastq", required=True, help="Input FASTQ file with UMI data")
    parser.add_argument("--seq_fastq", required=True, help="Input FASTQ file with sequences")
    parser.add_argument("--output_fastq", required=True, help="Output FASTQ file with cluster_id appended")
    args = parser.parse_args()

def extract_cluster_id(description):
    """Extracts the cluster_id from a FASTQ record description.
       Expects a field "cluster_id=XXXX".
    """
    for item in description.split():
        if item.startswith("cluster_id="):
            return item.split("=", 1)[1]
    raise ValueError(f"cluster_id not found in description: {description}")

def main():
    umi_iter = SeqIO.parse(args.umi_fastq, "fastq")
    seq_iter = SeqIO.parse(args.seq_fastq, "fastq")
    output_records = []
    
    for umi_rec, seq_rec in zip_longest(umi_iter, seq_iter):
        if umi_rec is None or seq_rec is None:
            sys.exit("Error: The two FASTQ files have different number of records.")
        try:
            cluster_id = extract_cluster_id(umi_rec.description)
        except ValueError as e:
            sys.exit(f"Error: {e}")
        # Append the cluster_id to the sequence record's description.
        seq_rec.description += f" unique_id={cluster_id}"
        output_records.append(seq_rec)
    
    SeqIO.write(output_records, args.output_fastq, "fastq")

if __name__ == '__main__':
    main()
