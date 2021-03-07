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
        timestampDir = lambda wildcards: config['timestamp']
    shell:
        """
        mkdir -p {params.timestampDir}
        if [ -d sequences/UMI ]; then
            find sequences/UMI -type f \(-name "*.csv" -o -name "*.tsv" \) | xargs mv -t {params.timestampDir}
            rm -r -f sequences/UMI
        fi
        find . -type f -wholename 'sequences/*.fastq.gz' -delete
        """

# clean up compute batches alignment
rule alignment_clean:
    input:
        [os.path.join(dirpath, f) for dirpath, _, files in os.walk('alignments') for f in files if f.endswith(('.bam','.bam.bai'))]
    output:
        touch('.alignment_clean.done')
    shell:
        "rm {input}"

# clean up compute batches demux
rule demux_clean:
    input:
        keep = [os.path.join(dirpath, f) for dirpath, _, files in os.walk('demux') for f in files if f.endswith('.csv')]
    output:
        touch('.demux_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp']
    shell:
        """
        if [ ! -z {input.keep} ]; then
            mkdir -p {params.timestampDir}
            mv {input.keep} -t {params.timestampDir}
        fi
        if [ -d demux ]; then
            rm -r demux
        fi
        """

rule mutation_data_clean:
    input:
        mutStats = [f for f in os.listdir('.') if f.endswith('_mutation-stats.csv')],
        keepDirs = [dirpath.strip('./') for dirpath, _, files in os.walk('.') if dirpath.endswith(('plots','mutSpectra','mutation_data'))]
    output:
        touch('.mutation_data_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp']
    shell:
        """
        if [ ! -z {input.keepDirs} ]; then
            mkdir -p {params.timestampDir}
            mv {input.keepDirs} -t {params.timestampDir}
        fi
        if [ ! -z {input.mutStats} ]; then
            mkdir -p {params.timestampDir}
            mv {input.mutStats} -t {params.timestampDir}
        fi
        """

rule logs_clean:
    output:
        touch('.logs_clean.done')
    params:
        timestampDir = lambda wildcards: config['timestamp']
    shell:
        """
        if [ -d log ]; then
            mkdir -p {params.timestampDir}/mapleLogs
            mv log/* {params.timestampDir}/mapleLogs
        fi
        if [ -d .snakemake/log ]; then
            mkdir -p {params.timestampDir}/snakemakeLogs
            mv .snakemake/log/* {params.timestampDir}/snakemakeLogs
        fi
        cp *.yaml {params.timestampDir}
        """

# clean up everything
rule clean:
    input:
        rules.sequences_clean.output,
        rules.alignment_clean.output,
        rules.demux_clean.output,
        rules.mutation_data_clean.output,
        rules.logs_clean.output
    shell:
        "rm {input}"

# Copyright (c) 2018-2020, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann, modified by Gordon Rix
# ---------------------------------------------------------------------------------