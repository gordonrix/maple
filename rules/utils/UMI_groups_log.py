#!/usr/bin/env python3
"""
UMI_groups_log.py

Parses a FASTQ file containing UMI data, extracts UMI‐related tags from each record’s description,
and builds a CSV log with the following columns:
    read_id, umi, umi_count, final_umi, final_umi_count, unique_id

It expects the FASTQ headers to include fields such as:
  - same_umi: the count of reads with the same exact UMI,
  - cluster_size: the count for the consensus (final) read of a cluster,
  - cluster_id: a unique identifier for the UMI group.

If any record is missing the "cluster_id" tag, the script exits with an error.
The final DataFrame is sorted by unique_id and written to the specified CSV file.
"""

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import argparse

# Argument parsing at the module level.
if 'snakemake' in globals():
    # When run via Snakemake, assume the input and output are provided in snakemake.input/output.
    class Args:
        pass
    args = Args()
    args.input_fastq = snakemake.input.fastq
    args.output_csv = snakemake.output.csv
else:
    parser = argparse.ArgumentParser(description="Generate UMI groups log CSV from FASTQ file")
    parser.add_argument("--input_fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output_csv", required=True, help="Output CSV file")
    args = parser.parse_args()

def description_to_dict(description):
    d = {}
    # Skip the first token (assumed to be the read name), then parse key=value pairs.
    for item in description.split(' ')[1:]:
        key, value = item.split('=')
        d[key] = value
    return d

def main():
    columns = ['read_id', 'umi', 'umi_count', 'final_umi_count', 'unique_id']
    rows = []
    count = 0

    for record in SeqIO.parse(args.input_fastq, "fastq"):
        desc_dict = description_to_dict(record.description)
        if 'cluster_id' not in desc_dict:
            sys.exit(f"Error: Record {record.id} is missing the cluster_id tag in its description.")
        umi = str(record.seq)  # Ensure the UMI is a plain string
        umi_count = desc_dict.get('same_umi', pd.NA)
        final_umi_count = desc_dict.get('cluster_size', pd.NA)
        unique_id = int(desc_dict['cluster_id'])
        count += 1
        rows.append([record.id, umi, umi_count, final_umi_count, unique_id])

    df = pd.DataFrame(rows, columns=columns)

    # helper columns for sorting
    df["is_na_fumi"] = df["final_umi_count"].isna().astype(int)
    df["is_na_umi"] = df["umi_count"].isna().astype(int)

    # Sort to enable consensus filling:
    # 1) unique_id ascending -> groups all rows of the same ID together
    # 2) is_na_fumi ascending -> non-NA (consensus) final_umi_count (0) appear before NA (1), but all consensus rows are tied
    # 3) final_umi_count descending -> bigger counts come before smaller
    # 4) umi ascending -> groups identical UMIs together alphabetically
    # 5) is_na_umi ascending -> non-NA (consensus) umi_count appear before NA, but all consensus rows are tied
    # 6) umi_count descending -> among non-NA values, bigger counts come first
    df = df.sort_values(
        by=["unique_id", "is_na_fumi", "final_umi_count", "umi", "is_na_umi", "umi_count"],
        ascending=[True, True, False, True, True, False],
        na_position="last",
        kind="stable"
    ).reset_index(drop=True)
    df = df.drop(columns=["is_na_fumi", "is_na_umi"])

    # Fill consensus values for final_umi_count and umi_count:
    df["final_umi_count_consensus"] = df["final_umi_count"].ffill().bfill()
    df["umi_count_consensus"] = df["umi_count"].ffill().bfill()

    # validate that non-NA values equal the consensus
    if not (df.loc[df["final_umi_count"].notna(), "final_umi_count"] == 
            df.loc[df["final_umi_count"].notna(), "final_umi_count_consensus"]).all():
        raise ValueError("Inconsistent final_umi_count values within a group")
    
    if not (df.loc[df["umi_count"].notna(), "umi_count"] ==
            df.loc[df["umi_count"].notna(), "umi_count_consensus"]).all():
        raise ValueError("Inconsistent umi_count values within a group")

    df["final_umi_count"] = df["final_umi_count_consensus"]
    df["umi_count"] = df["umi_count_consensus"]
    df.drop(columns=["final_umi_count_consensus", "umi_count_consensus"], inplace=True)
    df = df.sort_values(by=['final_umi_count', 'umi_count', 'unique_id'], ascending=[False, False, True])
    df.to_csv(args.output_csv, index=False)

if __name__ == '__main__':
    main()
