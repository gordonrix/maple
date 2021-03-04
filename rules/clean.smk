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
import os, sys, glob, datetime
# local rules
localrules: sequences_clean, alignment_clean, methylation_clean, sv_clean, demux_clean, transcript_isoforms_clean, clean

onstart:
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


# remove flag files to force re-run of cleanup
for f in glob.glob('.*.done'):
    os.remove(f)


# clean up compute batches basecalling, keep UMI stats files
rule sequences_clean:
    input:
        keep = [f for f in os.listdir(os.path.join('sequences', 'UMI')) if f.endswith(('.csv','.tsv'))],
        UMIdir = os.path.join('sequences', 'UMI')
        batches = [dirpath for dirpath, _, files in os.walk('sequences') if dirpath.endswith('batches')],
        consensus = [f for f in os.listdir('sequences') if f.endswith('_UMIconsensus.fastq.gz')]
    output:
        touch('.sequences_clean.done')
    params:
        timestampDir = lambda wildcards: timestamp
    shell:
        """
        mkdir -p {params.timestampDir}
        mv {input.keep} {params.timestampDir}/
        rm -r {input.UMIdir}
        rm -r {input.batches}
        rm {input.consensus)
        """

# clean up compute batches alignment
rule alignment_clean:
    input:
        directory('alignments')
    output:
        touch('.alignment_clean.done')
    shell:
        "rm -r {input}"

# clean up compute batches demux
rule demux_clean:
    input:
        keep = [f for f in os.listdir('demux') if f.endswith('.csv')],
        demux = directory('demux')
    output:
        touch('.demux_clean.done')
    params:
        timestampDir = lambda wildcards: timestamp
    shell:
        """
        mkdir -p {params.timestampDir}
        mv {input.keep} {params.timestampDir}/
        rm -r {input.dir}
        """

rule mutation_data_clean:
    input:
        keepDirs = 'plots mutSpectra mutation_data',
        mutStats = [f for f in os.listdir('.') if f.endswith('_mutation-stats.csv')]
    output:
        touch('.mutation_data_clean.done')
    params:
        timestampDir = lambda wildcards: timestamp
    shell:
        """
        mkdir -p {params.timestampDir}
        mv {input.keepDirs) {params.timestampDir}/
        mv {input.mutStats) {params.timestampDir}/
        """

# clean up everything
rule clean:
    input:
        rules.sequences_clean.output,
        rules.alignment_clean.output,
        rules.demux_clean.output,
        rules.mutation_data_clean.output
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