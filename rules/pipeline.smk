# imports
import os, sys, glob
from rules.utils.get_file import get_sequence_batch
from rules.utils.storage import get_flowcell, get_kit, get_ID
# local rules
if not config['merge_paired_end']:
    localrules: basecaller_merge_tag

# get batch of reads as IDs or fast5
def get_signal_batch(wildcards, config):
    raw_dir = config['storage_data_raw']
    batch_file = os.path.join(raw_dir, wildcards.runname, 'reads', wildcards.batch)
    if os.path.isfile(batch_file + '.tar'):
        return batch_file + '.tar'
    elif os.path.isfile(batch_file + '.fast5'):
        return batch_file + '.fast5'
    else:
        return []

# prefix of raw read batches
def get_batch_ids_raw(config, tag, runname):
    batches_tar, = glob_wildcards("{datadir}/{runname}/reads/{{id}}.tar".format(datadir=config["storage_data_raw"], runname=runname))
    batches_fast5, = glob_wildcards("{datadir}/{runname}/reads/{{id}}.fast5".format(datadir=config["storage_data_raw"], runname=runname))
    return batches_tar + batches_fast5

def get_batch_ids_sequences(config, tag, runname):
    batches_fastqgz, = glob_wildcards("sequences/batches/{tag}/{runname}/{{id}}.fastq.gz".format(tag=tag, runname=runname))
    return batches_fastqgz

# get batches
def get_batches_basecaller(wildcards):
    output = []
    for runname in config['runs'][wildcards.tag]['runname']:
        if config['do_basecalling']:
            outputs = expand("sequences/batches/{tag}/{runname}/{batch}.fastq.gz",
                                tag = wildcards.tag,
                                runname=runname,
                                batch=get_batch_ids_raw(config, wildcards.tag, runname))
        else:
            outputs = expand("sequences/batches/{tag}/{runname}/{batch}.fastq.gz",
                                tag = wildcards.tag,
                                runname = runname,
                                batch = get_batch_ids_sequences(config, tag, runname))
        output.extend(outputs)
    return output

if config['do_basecalling']:

    # guppy basecalling
    rule guppy:
        input:
            batch = lambda wildcards : get_signal_batch(wildcards, config),
            run = lambda wildcards : [os.path.join(config['storage_data_raw'], wildcards.runname)] + ([os.path.join(config['storage_data_raw'], wildcards.runname, 'reads.fofn')] if get_signal_batch(wildcards, config).endswith('.txt') else [])
        output:
            ["sequences/batches/{tag}/{runname}/{batch, [^.]*}.fastq.gz"] +
            ["sequences/batches/{tag}/{runname}/{batch, [^.]*}.sequencing_summary.txt"] +
            (["sequences/batches/{tag}/{runname}/{batch, [^.]*}.hdf5"] if config.get('basecalling_guppy_config') and 'modbases' in config['basecalling_guppy_config'] else [])
        shadow: "shallow"
        threads: config['threads_basecalling']
        resources:
            threads = lambda wildcards, threads: threads,
            mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.1 * (attempt - 1))) * (config['memory']['guppy_basecaller'][0] + config['memory']['guppy_basecaller'][1] * threads)),
            time_min = lambda wildcards, threads, attempt: int((1440 / threads) * attempt * config['runtime']['guppy_basecaller']), # 90 min / 16 threads
            gpu = 1
        params:
            guppy_config = lambda wildcards : '-c {cfg}{flags}'.format(
                                cfg = config.get('basecalling_guppy_config') or 'dna_r9.4.1_450bps_fast.cfg',
                                flags = ' --fast5_out' if config.get('basecalling_guppy_config') and 'modbases' in config['basecalling_guppy_config'] else ''),
            guppy_server = lambda wildcards, input : '' if (config.get('basecalling_guppy_flags') and '--port' in config['basecalling_guppy_flags']) else '--port ' + config['basecalling_guppy_server'][hash(input.batch) % len(config['basecalling_guppy_server'])] if config.get('basecalling_guppy_server') else '',
            guppy_flags = lambda wildcards : config.get('basecalling_guppy_flags') or '',
            filtering = lambda wildcards : '--qscore_filtering --min_qscore {score}'.format(score = config['basecalling_guppy_qscore_filter']) if config['basecalling_guppy_qscore_filter'] > 0 else '',
            index = lambda wildcards : '--index ' + os.path.join(config['storage_data_raw'], wildcards.runname, 'reads.fofn') if get_signal_batch(wildcards, config).endswith('.txt') else '',
            mod_table = lambda wildcards, input, output : output[2] if len(output) == 3 else ''
        # singularity:
        #     config['singularity_images']['basecalling']
        shell:
            """
            mkdir -p raw
            {config[bin_singularity][python]} {config[sbin_singularity][storage_fast5Index.py]} extract {input.batch} raw/ {params.index} --output_format bulk
            {config[bin_singularity][guppy_basecaller]} -i raw/ --recursive --num_callers 1 --cpu_threads_per_caller {threads} -s workspace/ {params.guppy_config}  {params.filtering} {params.guppy_flags} {params.guppy_server}
            FASTQ_DIR='workspace/pass'
            if [ \'{params.filtering}\' = '' ]; then
                FASTQ_DIR='workspace'
            fi
            find ${{FASTQ_DIR}} -regextype posix-extended -regex '^.*f(ast)?q' -exec cat {{}} \; | gzip > {output[0]}
            find ${{FASTQ_DIR}} -name 'sequencing_summary.txt' -exec mv {{}} {output[1]} \;
            """
    
rule basecaller_stats:
    input:
        lambda wildcards: get_batches_basecaller(wildcards)
    output:
        "sequences/batches/{tag}/{runname}.hdf5"
    run:
        import gzip
        import pandas as pd
        def fastq_iter(iterable):
            while True:
                try:
                    title = next(iterable)
                    assert title[0] == '@'
                    seq = next(iterable)
                    _ = next(iterable)
                    qual = next(iterable)
                except StopIteration:
                    return
                mean_q = sum([ord(x) - 33 for x in qual]) / len(qual) if qual else 0.0
                yield len(seq), mean_q
        line_iter = (line for f in input for line in gzip.open(f, 'rb').read().decode('utf-8').split('\n') if line)
        df = pd.DataFrame(fastq_iter(line_iter), columns=['length', 'quality'])
        df.to_hdf(output[0], 'stats')


if config['merge_paired_end']:
    rule merge_paired_end:
        input:
            fwd = lambda wildcards: os.path.join('sequences', 'paired', config['runs'][wildcards.tag]['fwdReads']),
            rvs = lambda wildcards: os.path.join('sequences', 'paired', config['runs'][wildcards.tag]['rvsReads'])
        output:
            protected("sequences/{tag, [^\/_]*}.fastq.gz")
        params:
            flags = config['NGmerge_flags']
        shell:
            """
            NGmerge -1 {input.fwd} -2 {input.rvs} -o {output} {params.flags}
            """

else:
    rule basecaller_merge_tag: # combine batches of basecalled reads into a single file
        input:
            lambda wildcards: get_batches_basecaller(wildcards)
        output:
            protected("sequences/{tag, [^\/_]*}.fastq.gz")
        run:
            with open(output[0], 'wb') as fp_out:
                for f in input:
                    with open(f, 'rb') as fp_in:
                        fp_out.write(fp_in.read())

rule UMI_minimap2:
    input:
        sequence = 'sequences/{tag}.fastq.gz',
        reference = lambda wildcards: config['runs'][wildcards.tag]['reference']
    output:
        pipe("sequences/UMI/{tag, [^\/_]*}_noConsensus.sam")
    threads: config['threads_alignment']
    group: "minimap2"
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    params:
        tempRef = lambda wildcards: wildcards.tag.split('.f')[0] + 'tempRef.fasta'
    # singularity:
    #     config['singularity_images']['alignment']
    shell:
        """
        head -2 {input.reference} > {params.tempRef}
        {config[bin_singularity][minimap2]} -t {threads} {config[alignment_minimap2_flags]} {params.tempRef} {input.sequence} 1>> {output} 2> >(tee {output}.log >&2)
        if [ $(grep 'ERROR' {output}.log | wc -l) -gt 0 ]; then exit 1; else rm {output}.log; fi
        rm {params.tempRef}
        """

# sam to bam conversion
rule UMI_aligner_sam2bam:
    input:
        sam = "sequences/UMI/{tag}_noConsensus.sam"
    output:
        bam = temp("sequences/UMI/{tag, [^\/_]*}_noConsensus.bam"),
        index = temp("sequences/UMI/{tag, [^\/_]*}_noConsensus.bam.bai")
    shadow: "minimal"
    threads: 1
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: int((1.0 + (0.2 * (attempt - 1))) * 5000)
    # singularity:
    #     config['singularity_images']['alignment']
    shell:
        """
        {config[bin_singularity][samtools]} view -b {input.sam} | {config[bin_singularity][samtools]} sort -m 4G > {output.bam}
        {config[bin_singularity][samtools]} index {output.bam}
        """

rule UMI_extract:
    input:
        bam = 'sequences/UMI/{tag}_noConsensus.bam',
        index = "sequences/UMI/{tag}_noConsensus.bam.bai"
    output:
        extracted = temp('sequences/UMI/{tag, [^\/_]*}_UMIextract.bam'),
        index = temp('sequences/UMI/{tag, [^\/_]*}_UMIextract.bam.bai'),
        log = 'sequences/UMI/{tag, [^\/_]*}_UMI-extract.csv'
    params:
        barcode_contexts = lambda wildcards: [config['runs'][wildcards.tag]['barcodeInfo'][barcodeType]['context'].upper() for barcodeType in config['runs'][wildcards.tag]['barcodeInfo']] if config['demux'] else None,
        reference = lambda wildcards: config['runs'][wildcards.tag]['reference'],
        UMI_contexts = lambda wildcards: config['runs'][wildcards.tag]['UMI_contexts']
    script:
        'utils/UMI_extract.py'

rule UMI_group:
    input:
        bam = 'sequences/UMI/{tag}_UMIextract.bam',
        index = 'sequences/UMI/{tag}_UMIextract.bam.bai'
    output:
        bam = temp('sequences/UMI/{tag, [^\/_]*}_UMIgroup.bam'),
        index = temp('sequences/UMI/{tag, [^\/_]*}_UMIgroup.bam.bai'),
        log = 'sequences/UMI/{tag, [^\/_]*}_UMIgroup-log.tsv'
    shell:
        """
        umi_tools group -I {input.bam} --group-out={output.log} --output-bam --per-gene --gene-tag=GN --edit-distance-threshold 2 -S {output.bam}
        samtools index {output.bam}
        """

rule plot_UMI_group:
    input:
        'sequences/UMI/{tag}_UMIgroup-log.tsv'
    output:
        csv = 'sequences/UMI/{tag, [^\/_]*}_UMIgroup-distribution.csv',
        plot = 'plots/{tag, [^\/_]*}_UMIgroup-distribution.html'
    script:
        'utils/plot_UMI_groups_distribution.py'

rule UMI_consensus:
    input:
        grouped = 'sequences/UMI/{tag}_UMIgroup.bam',
        index = 'sequences/UMI/{tag}_UMIgroup.bam.bai',
        log = 'sequences/UMI/{tag}_UMIgroup-log.tsv'
    output:
        temp('sequences/UMI/{tag, [^\/_]*}_UMIconsensus.fastq')
    script:
        'utils/UMI_consensus.py'

rule UMI_compress:
    input:
        'sequences/UMI/{tag}_UMIconsensus.fastq'
    output:
        'sequences/{tag, [^\/_]*}_UMIconsensus.fastq.gz'
    params:
        temp = lambda wildcards: f'sequences/UMI/{wildcards.tag}_UMIconsensus.fastq.gz'
    shell:
        """
        gzip {input}
        mv {params.temp} {output}
        """

def alignment_sequence_input(wildcards):
    if config['UMI_consensus']:
        return 'sequences/{tag}_UMIconsensus.fastq.gz'
    else:
        return 'sequences/{tag}.fastq.gz'

rule minimap2:
    input:
        sequence = alignment_sequence_input,
        reference = lambda wildcards: config['runs'][wildcards.tag]['reference']
    output:
        pipe("alignments/{tag, [^\/_]*}.sam")
    threads: config['threads_alignment']
    group: "minimap2"
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    params:
        tempRef = lambda wildcards: wildcards.tag.split('.f')[0] + 'tempRef.fasta'
    # singularity:
    #     config['singularity_images']['alignment']
    shell:
        """
        head -2 {input.reference} > {params.tempRef}
        {config[bin_singularity][minimap2]} -t {threads} {config[alignment_minimap2_flags]} {params.tempRef} {input.sequence} 1>> {output} 2> >(tee {output}.log >&2)
        if [ $(grep 'ERROR' {output}.log | wc -l) -gt 0 ]; then exit 1; else rm {output}.log; fi
        rm {params.tempRef}
        """

# sam to bam conversion
rule aligner_sam2bam:
    input:
        sam = "alignments/{tag}.sam"
    output:
        bam = "alignments/{tag, [^\/_]*}.bam",
        bai = "alignments/{tag, [^\/_]*}.bam.bai"
    shadow: "minimal"
    threads: 1
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: int((1.0 + (0.2 * (attempt - 1))) * 5000)
    # singularity:
    #     config['singularity_images']['alignment']
    shell:
        """
        {config[bin_singularity][samtools]} view -b {input.sam} | {config[bin_singularity][samtools]} sort -m 4G > {output.bam}
        {config[bin_singularity][samtools]} index {output.bam}
        """

# mapping stats
rule aligner_stats:
    input:
        "alignments/batches/{tag}/{runname}_aligned-file-names.txt.txt"
    output:
        "alignments/batches/{tag, [^\/]*}/{runname, [^.\/]*}.hdf5"
    threads: config.get('threads_samtools') or 1
    resources:
        threads = lambda wildcards, threads: threads,
    # singularity:
    #     config['singularity_images']['alignment']
    shell:
        """
        while IFS= read -r bam_file; do {config[bin_singularity][samtools]} view ${{bam_file}}; done < {input} | {config[bin_singularity][python]} {config[sbin_singularity][alignment_stats.py]} {output}
        """

if config['demux']:

    rule generate_barcode_ref:
        input:
            'alignments/{tag}.bam'
        output:
            temp('demux/{tag, [^\/_]*}_generate_barcode_ref_complete.txt')
        script:
            'utils/generate_barcode_ref.py'

    checkpoint demultiplex:
        input:
            aln = 'alignments/{tag}.bam',
            flag = 'demux/{tag}_generate_barcode_ref_complete.txt'
        output:
            temp('demux/{tag, [^\/_]*}_demultiplex_complete.txt')
        script:
            'utils/demux.py'

    rule index_demuxed:
        input:
            'demux/{tag}_{barcodes}.bam'
        output:
            'demux/{tag, [^\/]*}_{barcodes}.bam.bai'
        shell:
            """
            samtools index {input}
            """

    rule mutation_analysis:
        input:
            bam = 'demux/{tag}_{barcodes}.bam',
            bai = 'demux/{tag}_{barcodes}.bam.bai'
        output:
            expand('mutation_data/{{tag, [^\/]*}}_{{barcodes}}_{datatype}', datatype = ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-aggregated.csv', 'AA-muts-distribution.csv'] if config['do_AA_analysis'] else ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv'])
        script:
            'utils/mutation_analysis.py'

else:
    rule mutation_analysis:
        input:
            bam = 'alignments/{tag}.bam',
            bai = 'alignments/{tag}.bam.bai'
        output:
            expand('mutation_data/{{tag, [^\/]*}}_{{barcodes}}_{datatype}', datatype = ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-aggregated.csv', 'AA-muts-distribution.csv'] if config['do_AA_analysis'] else ['highest-abundance-alignments.txt', 'genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv'])
        script:
            'utils/mutation_analysis.py'

def mut_stats_input(wildcards):
    datatypes = ['genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv', 'AA-muts-aggregated.csv', 'AA-muts-distribution.csv'] if config['do_AA_analysis'] else ['genotypes.csv', 'failures.csv', 'NT-muts-aggregated.csv', 'NT-muts-distribution.csv']
    if config['demux']:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex_complete')[0]
        checkpoint_demux_files = checkpoint_demux_prefix + '{BCs}.bam'
        return expand('mutation_data/{tag}_{barcodes}_{datatype}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, datatype=datatypes)
    else:
        return expand('mutation_data/{tag}_all_{datatype}', tag=wildcards.tag, datatype=datatypes)

rule mut_stats:
	input:
		mut_stats_input
	output:
		'{tag, [^\/]*}_mutation-stats.csv'
	script:
		'utils/mutation_statistics.py'

rule plot_mutation_spectrum:
    input:
        '{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/]*}_mutation-spectra.html'
    script:
        'utils/plot_mutation_spectrum.py'

def plot_mutations_aggregated_input(wildcards):
    if config['demux']:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex_complete')[0]
        checkpoint_demux_files = checkpoint_demux_prefix + '{BCs}.bam'
        return expand('mutation_data/{tag}_{barcodes}_{AAorNT}-muts-aggregated.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, AAorNT = wildcards.AAorNT)
    else:
        return 'mutation_data/{tag}_all_{AAorNT}-muts-aggregated.csv'

rule plot_mutations_aggregated:
    input:
        aggregated = plot_mutations_aggregated_input,
        mutStats = '{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/]*}_{AAorNT}-mutations-aggregated.html'
    script:
        'utils/plot_mutations_aggregated.py'

def plot_mutations_distribution_input(wildcards):
    if config['demux']:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex_complete')[0]
        checkpoint_demux_files = checkpoint_demux_prefix + '{BCs}.bam'
        return expand('mutation_data/{tag}_{barcodes}_{AAorNT}-muts-distribution.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, AAorNT = wildcards.AAorNT)
    else:
        return 'mutation_data/{tag}_all_{AAorNT}-muts-distribution.csv'

rule plot_mutations_distribution:
    input:
        plot_mutations_distribution_input
    output:
        'plots/{tag, [^\/]*}_{AAorNT}-mutation-distributions.html'
    script:
        'utils/plot_mutation_distribution.py'


#         """

# ---------------------------------------------------------------------------------
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