# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : remove large files and move analysis files to timestamped folder
#                   so that analysis may be repeated
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------

# imports
import os, sys, glob
from datetime import datetime
# local rules
localrules: sequences_clean, alignment_clean, demux_clean, mutation_data_clean, logs_clean, clean

config['timestamp'] = datetime.now().strftime("%Y%m%d-%H%M%S")

# remove flag files to force re-run of cleanup
for f in glob.glob('.*.done'):
    os.remove(f)


# clean up sequences folder: remove combined .fastq.gz files (but not basecalled batches), move UMI stats files
rule sequences_clean:
    output:
        touch('.sequences_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp'] + '-mapleOutputs'
    shell:
        """
        mkdir -p {params.timestamp_dir}
        if [ -d sequences ]; then
            if [ -d sequences/UMI ]; then
                find sequences/UMI -type f \( -name "UMI-extract-summary.csv" \) -exec cp {{}} {params.timestamp_dir} \;
                find sequences/UMI -type f \( -name "*.fasta.gz" \) -exec cp {{}} {params.timestamp_dir} \;
            fi
            if [ -d sequences/paired ]; then
                find sequences/paired -type f \( -name "*_failed-merge_1.fastq.gz" \) -delete
                find sequences/paired -type f \( -name "*_failed-merge_2.fastq.gz" \) -delete
                find sequences/paired -type f \( -name "*_NGmerge.log" \) -delete
            fi
        fi
        """

# clean up compute batches alignment
rule alignment_clean:
    output:
        touch('.alignment_clean.done')
    shell:
        """
        if [ -d alignments ]; then
            rm -r alignments
        fi
        """

# clean up compute batches demux
rule demux_clean:
    input:
        '.mutation_data_clean.done'
    output:
        touch('.demux_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp'] + '-mapleOutputs'
    shell:
        """
        if [ -d demux ]; then
            rm -r demux
        fi
        """

rule mutation_data_clean:
    output:
        touch('.mutation_data_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp'] + '-mapleOutputs',
        keep = [directoryORfile for directoryORfile in os.listdir('.') if directoryORfile in ['plots', 'mutSpectra', 'mutation_data', 'mutation-stats.csv', 'demux-stats.csv', 'dms-view-table.csv', 'maple']]
    shell:
        """
        if [ ! -z "{params.keep}" ]; then
            mkdir -p {params.timestamp_dir}
            mv {params.keep} -t {params.timestamp_dir}
        fi
        """

rule logs_clean:
    output:
        touch('.logs_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp'] + '-mapleOutputs',
        metadata_dir = lambda wildcards: config['metadata']
    shell:
        """
        if [ -d log ]; then
            mkdir -p {params.timestamp_dir}/mapleLogs
            mv log/* {params.timestamp_dir}/mapleLogs
        fi
        if [ -d .snakemake/log ]; then
            mkdir -p {params.timestamp_dir}/snakemakeLogs
            mv .snakemake/log/* {params.timestamp_dir}/snakemakeLogs
        fi
        cp *.yaml {params.timestamp_dir}
        cp -r {params.metadata_dir} {params.timestamp_dir}
        """

# clean up everything
rule clean:
    input:
        rules.sequences_clean.output,
        rules.alignment_clean.output,
        rules.demux_clean.output,
        rules.mutation_data_clean.output,
        rules.logs_clean.output
    params:
        timestamp_dir = lambda wildcards: config['timestamp'] + '-mapleOutputs'
    shell:
        """
        rm {input}
        zip -r -m {params.timestamp_dir}.zip {params.timestamp_dir}
        """
