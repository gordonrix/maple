#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Rules for retrieval of sequence files from locations outside of the working directory
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

import os
import glob
import pandas as pd
import subprocess
import re

from rules.utils.common import retrieve_fastqs
    

rule combine_tag: # combine batches of reads into a single file
    input:
        lambda wildcards: retrieve_fastqs(config['sequences_dir'], config['runs'][wildcards.tag].get('runname',''), config['fastq_dir'])
    output:
        temp("sequences/{tag, [^\/_]*}_combined.fastq.gz")
    run:
        with open(output[0], 'wb') as fp_out:
            if len(input)==0:
                raise RuntimeError(f"Basecalled sequence batches not found for tag `{wildcards.tag}`.")
            for f in input:
                with open(f, 'rb') as fp_in:
                    fp_out.write(fp_in.read())

checkpoint basespace_retrieve_project:
    output:
        flag = "sequences/bs{bs_project_ID}/.done"
    shell:
        """
        bs download project -i {wildcards.bs_project_ID} -o sequences/bs{wildcards.bs_project_ID} --extension fastq.gz
        touch {output.flag}
        """

def get_bs_paired_end(tag):
    """
    Retrieve paired-end fastq files from BaseSpace for a given tag.
    Assumes paired end files are located in the following directory structure after a bs project download:
        sequences/bs{bs_project_ID}/{sample_ID}*/{tag}*R1*.fastq.gz
        sequences/bs{bs_project_ID}/{sample_ID}*/{tag}*R2*.fastq.gz
    """
    project_ID = config['runs'][tag].get('bs_project_ID', None)
    if not project_ID:
        raise ValueError(f"No BaseSpace project ID found for tag {tag}")
    flag = checkpoints.basespace_retrieve_project.get(bs_project_ID=project_ID).output[0]
    sample_ID = config['runs'][tag].get('sample_ID', None)
    if not sample_ID:
        raise ValueError(f"[ERROR] Sample ID not provided for tag {tag}, but bs_project_ID is provided.")
    project_dir = os.path.dirname(flag)
    fwd = glob.glob(os.path.join(project_dir, f"{config['runs'][tag]['sample_ID']}*", f"*R1*.fastq.gz"))
    rvs = glob.glob(os.path.join(project_dir, f"{config['runs'][tag]['sample_ID']}*", f"*R2*.fastq.gz"))
    if fwd and rvs:
        return fwd[0], rvs[0]
    else:
        raise FileNotFoundError(f"Paired-end files for sample {tag} not found in {project_dir}")

rule merge_paired_end:
    input:
        fwd = lambda wildcards: retrieve_fastqs(config['sequences_dir'], [config['runs'][wildcards.tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][wildcards.tag]['fwdReads'])[0],
        rvs = lambda wildcards: retrieve_fastqs(config['sequences_dir'], [config['runs'][wildcards.tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][wildcards.tag]['rvsReads'])[0]
    output:
        merged = "sequences/paired/{tag, [^\/_]*}.fastq.gz",
        log = "sequences/paired/{tag, [^\/_]*}_NGmerge.log",
        failedfwd = "sequences/paired/{tag, [^\/_]*}_failed-merge_1.fastq.gz",
        failedrvs = "sequences/paired/{tag, [^\/_]*}_failed-merge_2.fastq.gz"
    params:
        flags = config['NGmerge_flags'],
        failed = "sequences/paired/{tag, [^\/_]*}_failed-merge"
    shell:
        """
        NGmerge -1 {input.fwd} -2 {input.rvs} -o {output.merged} -l {output.log} -f {params.failed} -z {params.flags}
        """

rule merge_paired_end_from_bs:
    input:
        fwd = lambda wildcards: get_bs_paired_end(wildcards.tag)[0],
        rvs = lambda wildcards: get_bs_paired_end(wildcards.tag)[1]
    output:
        merged = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}.fastq.gz",
        log = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_NGmerge.log",
        failedfwd = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_failed-merge_1.fastq.gz",
        failedrvs = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_failed-merge_2.fastq.gz"
    params:
        flags = config['NGmerge_flags'],
        failed = "sequences/bs{bs_project_ID, [^\/_]*}_merged/{tag, [^\/_]*}_failed-merge"
    shell:
        """
        NGmerge -1 {input.fwd} -2 {input.rvs} -o {output.merged} -l {output.log} -f {params.failed} -z {params.flags}
        """

def get_merged_seqs(tag):

    if any([x in config['runs'][tag] for x in ['fwdReads', 'rvsReads']]):
        return f'sequences/paired/{tag}.fastq.gz'

    elif 'runname' in config['runs'][tag]:
        return f'sequences/{tag}_combined.fastq.gz'

    elif 'bs_project_ID' in config['runs'][tag]:
        project_ID = config['runs'][tag]['bs_project_ID']
        return f'sequences/bs{project_ID}_merged/{tag}.fastq.gz'
    
    else:
        return ''

rule move_seqs: # allows for merging batches of sequences or merging paired end reads depending on the tag definition
    input:
        lambda wildcards: get_merged_seqs(wildcards.tag)
    output:
        'sequences/{tag, [^\/_]*}.fastq.gz'
    shell:
        """
        mv {input} {output}
        """