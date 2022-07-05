# imports
import os, sys, glob

# local rules
if not config['merge_paired_end']:
    localrules: basecaller_merge_tag

# get batch of reads fast5
def get_signal_batch(wildcards):
    batch_file = os.path.join(config['minknowDir'], wildcards.expt, wildcards.sample, wildcards.runname, config['fast5_dir'], wildcards.batch + '.fast5')
    if os.path.isfile(batch_file):
        return batch_file
    else:
        print(f"[WARNING] Batch file `{batch_file}` was requested but not found.")
        return []

# prefix of raw read batches
def get_batch_ids_raw_minknow(runpath, fast5_dir):
    batches_fast5, = glob_wildcards("{runpath}/{fast5_dir}/{{id}}.fast5".format(runpath=runpath, fast5_dir=fast5_dir))
    return batches_fast5

def get_batch_ids_sequences(runname):
    batches_fastqgz, = glob_wildcards("sequences/batches/{runname}/{{id}}.fastq.gz".format(runname=runname))
    return batches_fastqgz

def get_batch_ids_sequences_minkow(runpath, fastq_dir):
    batches_fastqgz, = glob_wildcards("{runpath}/{fastq_dir}/{{id}}.fastq.gz".format(runpath=runpath, fastq_dir=fastq_dir))
    return batches_fastqgz

# get batches
def get_batches_basecaller(wildcards):
    output = []
    for runname in config['runs'][wildcards.tag]['runname']:
        outputs = []
        if not config['do_basecalling']: # basecalling was already performed so basecalled sequences are fetched
            outputs = expand("sequences/batches/{runname}/{batch}.fastq.gz",
                                runname = runname,
                                batch = get_batch_ids_sequences(runname))

        if outputs==[] and os.path.isdir(config['minknowDir']): # no sequence batches were found within the working directory, so will try to find them in the minknow directory instead
            experimentDirs = [d for d in os.listdir(config['minknowDir']) if os.path.isdir(os.path.join(config['minknowDir'],d))]
            expt_sample = False
            for expt in experimentDirs:
                if expt_sample: break
                sampleDirs = [d for d in os.listdir(os.path.join(config['minknowDir'], expt)) if os.path.isdir(os.path.join(config['minknowDir'], expt,d))]
                for sample in sampleDirs:
                    if runname in os.listdir(os.path.join(config['minknowDir'], expt, sample)):
                        expt_sample = (expt, sample)
                        break

            if expt_sample:
                if config['do_basecalling'] == True: # choose whether to request fastq files for all fast5 files or just what fastq files are already present based on whether basecalling should be performed
                    batch = get_batch_ids_raw_minknow(os.path.join(config['minknowDir'].rstrip('/'), expt_sample[0], expt_sample[1], runname), config['fast5_dir'])
                    outputs.extend(expand("sequences/batches/{expt}/{sample}/{runname}/{fastq_dir}/{batch}.fastq.gz",
                                            minknowDir = config['minknowDir'].rstrip('/'),
                                            expt = expt_sample[0],
                                            sample = expt_sample[1],
                                            runname = runname,
                                            fastq_dir = config['fastq_dir'],
                                            batch = batch))
                else:
                    batch = get_batch_ids_sequences_minkow(os.path.join(config['minknowDir'].rstrip('/'), expt_sample[0], expt_sample[1], runname), config['fastq_dir'])
                    outputs.extend(expand("{minknowDir}/{expt}/{sample}/{runname}/{fastq_dir}/{batch}.fastq.gz",
                                            minknowDir = config['minknowDir'].rstrip('/'),
                                            expt = expt_sample[0],
                                            sample = expt_sample[1],
                                            runname = runname,
                                            fastq_dir = config['fastq_dir'],
                                            batch = batch))
        output.extend(outputs)
    return output


if config['do_basecalling']:

    # guppy basecalling
    rule guppy:
        input:
            batch = lambda wildcards : get_signal_batch(wildcards)
        output:
            "sequences/batches/{expt}/{sample}/{runname}/{fastq_dir}/{batch, [^.]*}.fastq.gz"
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
            filtering = lambda wildcards : '--min_qscore {score}'.format(score = config['basecalling_guppy_qscore_filter']) if config['basecalling_guppy_qscore_filter'] > 0 else 0,
            mod_table = lambda wildcards, input, output : output[2] if len(output) == 3 else ''
        # singularity:
        #     config['singularity_images']['basecalling']
        shell:
            """
            mkdir -p raw
            {config[bin_singularity][python]} {config[sbin_singularity][storage_fast5Index.py]} extract {input.batch} raw/ --output_format bulk
            {config[bin_singularity][guppy_basecaller]} -i raw/ --recursive --num_callers 1 -s workspace/ {params.guppy_config} {params.filtering} {params.guppy_flags} {params.guppy_server}
            FASTQ_DIR='workspace/pass'
            if [ \'{params.filtering}\' = '' ]; then
                FASTQ_DIR='workspace'
            fi
            find ${{FASTQ_DIR}} -regextype posix-extended -regex '^.*f(ast)?q' -exec cat {{}} \; | gzip > {output[0]}
            """
    
rule basecaller_stats:
    input:
        lambda wildcards: get_batches_basecaller(wildcards)
    output:
        "sequences/batches/{runname}.hdf5"
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
            merged = "sequences/{tag, [^\/_]*}.fastq.gz",
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

elif config['nanopore']:
    rule basecaller_merge_tag: # combine batches of basecalled reads into a single file
        input:
            lambda wildcards: get_batches_basecaller(wildcards)
        output:
            "sequences/{tag, [^\/_]*}.fastq.gz"
        run:
            with open(output[0], 'wb') as fp_out:
                if len(input)==0:
                    raise RuntimeError(f"Basecalled sequence batches not found for tag `{wildcards.tag}`.")
                for f in input:
                    with open(f, 'rb') as fp_in:
                        fp_out.write(fp_in.read())

rule RCA_consensus:
    input:
        sequence = 'sequences/{tag}.fastq.gz',
        splintRef = lambda wildcards: os.path.join(config['references_directory'], f".{wildcards.tag}_splint.fasta")
    output:
        consensus = 'sequences/RCA-consensus/{tag, [^\/_]*}_RCA-consensus.fastq.gz',
        log = 'sequences/{tag, [^\/_]*}_RCA-consensus.log'
    threads: workflow.cores/len([tag for tag in config['do_RCA_consensus'] if config['do_RCA_consensus'][tag]==True])
    resources:
        threads = lambda wildcards, threads: threads,
    shell:
        """
        rm -r -f sequences/tempOutDir-{wildcards.tag}_RCA-consensus
        python3 -m C3POa -r {input.sequence} -o sequences/tempOutDir-{wildcards.tag}_RCA-consensus -s {input.splintRef} -n {threads} -co -z
        mv sequences/tempOutDir-{wildcards.tag}_RCA-consensus/splint/R2C2_Consensus.fasta.gz {output.consensus}
        mv sequences/tempOutDir-{wildcards.tag}_RCA-consensus/c3poa.log {output.log}
        rm -r sequences/tempOutDir-{wildcards.tag}_RCA-consensus
        
        """

rule UMI_minimap2:
    input:
        sequence = lambda wildcards: 'sequences/RCA-consensus/{tag, [^\/_]*}_RCA-consensus.fastq.gz' if config['do_RCA_consensus'][wildcards.tag] else 'sequences/{tag}.fastq.gz',
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        aln = pipe("sequences/UMI/{tag, [^\/_]*}_noConsensus.sam"),
        log = "sequences/UMI/{tag, [^\/_]*}_noConsensus.log"
    threads: config['threads_alignment']
    group: "minimap2"
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    shell:
        """
        {config[bin_singularity][minimap2]} -t {threads} {config[alignment_minimap2_flags]} {input.alnRef} {input.sequence} 1>> {output.aln} 2> >(tee {output.log} >&2)
        if [ $(grep 'ERROR' {output.log} | wc -l) -gt 0 ]; then exit 1; fi
        """

# sam to bam conversion
rule UMI_aligner_sam2bam:
    input:
        sam = "sequences/UMI/{tag}_noConsensus.sam"
    output:
        bam = "sequences/UMI/{tag, [^\/_]*}_noConsensus.bam",
        index = "sequences/UMI/{tag, [^\/_]*}_noConsensus.bam.bai"
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
        barcode_contexts = lambda wildcards: [config['runs'][wildcards.tag]['barcodeInfo'][barcodeType]['context'].upper() for barcodeType in config['runs'][wildcards.tag]['barcodeInfo']] if config['do_demux'][wildcards.tag] else None,
        reference = lambda wildcards: config['runs'][wildcards.tag]['reference'],
        UMI_contexts = lambda wildcards: config['runs'][wildcards.tag]['UMI_contexts']
    script:
        'utils/UMI_extract.py'

# collapse UMI extract log files into a single small file that summarizes UMI recognition in aggregate instead of on a read-by-read basis

rule UMI_extract_summary:
    input:
        expand('sequences/UMI/{tag}_UMI-extract.csv', tag = list(set( [config['consensusCopyDict'][str(t)] for t in config['runs'] if config['do_UMI_analysis'][t]] )))
    output:
        'sequences/UMI/UMI-extract-summary.csv'
    run:
        import pandas as pd
        maxUMIs = 0
        for tag in config['runs']:
            if 'UMI_contexts' in config['runs'][tag]:
                numUMIs = len(config['runs'][tag]['UMI_contexts'])
                if numUMIs > maxUMIs:
                    maxUMIs = numUMIs
        UMIcols = [f'umi_{UMInum}_failure' for UMInum in range(1,maxUMIs+1)]
        outDF = pd.DataFrame(columns=['tag','success','failure']+UMIcols)
        for f in input:
            df = pd.read_csv(f, dtype={'umi':str})
            dfCols = df.columns[2:]
            tag = f.split('/')[-1].split('_')[0]
            df['tag'] = tag
            dfSum = df.groupby(['tag'])[dfCols].sum().reset_index()
            if len(outDF.columns) != len(dfSum.columns): # add columns to match the maximum number of UMIs in a tag to allow for df concatenation
                for col in UMIcols:
                    if col not in dfCols:
                        dfSum[col] = 'NA'
            outDF = pd.concat([outDF, dfSum], ignore_index=True)
        outDF.to_csv(output[0])

rule UMI_group:
    input:
        bam = 'sequences/UMI/{tag}_UMIextract.bam',
        index = 'sequences/UMI/{tag}_UMIextract.bam.bai'
    output:
        bam = 'sequences/UMI/{tag, [^\/_]*}_UMIgroup.bam',
        index = 'sequences/UMI/{tag, [^\/_]*}_UMIgroup.bam.bai',
        log = 'sequences/UMI/{tag, [^\/_]*}_UMIgroup-log.tsv'
    params:
        UMI_mismatches = lambda wildcards: config['UMI_mismatches']
    shell:
        """
        umi_tools group -I {input.bam} --group-out={output.log} --output-bam --per-gene --gene-tag=GN --edit-distance-threshold {params.UMI_mismatches} -S {output.bam}
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

checkpoint UMI_splitBAMs:
    input:
        grouped = 'sequences/UMI/{tag}_UMIgroup.bam',
        index = 'sequences/UMI/{tag}_UMIgroup.bam.bai',
        log = 'sequences/UMI/{tag}_UMIgroup-log.tsv'
    output:
        outDir = directory('sequences/UMI/{tag, [^\/_]*}_UMIsplit')
    script:
        'utils/UMI_splitBAMs.py'

rule UMI_consensus:
    input:
        bam = 'sequences/UMI/{tag}_UMIsplit/UMI_{UMIID_finalUMI}_reads.bam',
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        hdf = 'sequences/UMI/{tag}_UMIsplit/UMI_{UMIID_finalUMI}_medaka_calls.hdf',
        consensus = 'sequences/UMI/{tag}_UMIsplit/UMI_{UMIID_finalUMI}_medaka_consensus.fasta'
    params:
        medaka_model = lambda wildcards: config['medaka_model'],
        refLength = lambda wildcards, input: len(next(SeqIO.parse(input.alnRef, 'fasta')).seq),
        consensus_flags = lambda wildcards: config['medaka_conensus_flags'],
        stitch_flags = lambda wildcards: config['medaka_stitch_flags']
    threads: lambda wildcards: config['threads_medaka']
    shell:
        """
        samtools sort -o {input.bam} {input.bam}
        samtools index {input.bam}
        medaka consensus {params.consensus_flags} --threads {threads} --chunk_len {params.refLength} --model {params.medaka_model} {input.bam} {output.hdf}
        medaka stitch {params.stitch_flags} --threads {threads} {output.hdf} {input.alnRef} {output.consensus}
        """

def UMI_merge_consensus_seqs_input(wildcards):
    UMI_split_output = checkpoints.UMI_splitBAMs.get(tag=wildcards.tag).output[0]
    UMI_split_output_template = os.path.join(UMI_split_output, 'UMI_{UMIID_finalUMI}_reads.bam')
    return expand('sequences/UMI/{tag}_UMIsplit/UMI_{UMIID_finalUMI}_medaka_consensus.fasta', tag=wildcards.tag, UMIID_finalUMI=glob_wildcards(UMI_split_output_template).UMIID_finalUMI)

rule UMI_merge_consensus_seqs:
    input:
        UMI_merge_consensus_seqs_input
    output:
        'sequences/UMI/{tag, [^\/_]*}_UMIconsensus.fasta'
    run:
        with open(output[0], 'w') as fp_out:
            if len(input)==0:
                raise RuntimeError(f"Consensus sequences not found for tag `{wildcards.tag}`.")
            for f in input:
                with open(f, 'r') as fp_in:
                    readID = f.split('/')[-1].strip('_medaka_consensus.fasta')
                    fp_out.write(f'>{readID}\n')
                    fp_out.write(fp_in.readlines()[1])

rule UMI_compress_consensus:
    input:
        'sequences/UMI/{tag}_UMIconsensus.fasta'
    output:
        'sequences/UMI/{tag, [^\/_]*}_UMIconsensus.fasta.gz'
    shell:
        """
        gzip {input}
        """

def alignment_sequence_input(wildcards):
    if config['do_UMI_analysis'][wildcards.tag]:
        return expand('sequences/UMI/{tag}_UMIconsensus.fasta.gz', tag=config['consensusCopyDict'][wildcards.tag])
    elif config['do_RCA_consensus'][wildcards.tag]:
        return 'sequences/RCA-consensus/{tag}_RCA-consensus.fastq.gz'
    else:
        return 'sequences/{tag}.fastq.gz'

rule minimap2:
    input:
        sequence = alignment_sequence_input,
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        aln = pipe("alignments/{tag, [^\/_]*}.sam"),
        log = "alignments/{tag, [^\/_]*}.log"
    threads: config['threads_alignment']
    group: "minimap2"
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    shell:
        """
        {config[bin_singularity][minimap2]} -t {threads} {config[alignment_minimap2_flags]} {input.alnRef} {input.sequence} 1>> {output.aln} 2> >(tee {output.log} >&2)
        if [ $(grep 'ERROR' {output.log} | wc -l) -gt 0 ]; then exit 1; fi
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

# Use NanoPlot package to generate and organize plots that describe sequencing data and alignments

rule NanoPlot_fastq:
    input:
        'sequences/{tag}.fastq.gz'
    output:
        'plots/nanoplot/{tag, [^\/_]*}_fastq_NanoStats.txt'
    params:
        fastqType = lambda wildcards: '--fastq_rich' if config['nanopore']==True else '--fastq',
        flags = lambda wildcards: config['nanoplot_flags']
    shell:
        """
        mkdir plots/nanoplot -p
        NanoPlot {params.flags} -o plots/nanoplot -p {wildcards.tag}_fastq_ {params.fastqType} {input[0]}
        """

rule NanoPlot_alignment_preConsensus:
    input:
        'sequences/UMI/{tag}_noConsensus.bam'
    output:
        'plots/nanoplot/{tag, [^\/_]*}_alignment_preConsensus_NanoStats.txt'
    params:
        flags = lambda wildcards: config['nanoplot_flags']
    shell:
        """
        mkdir plots/nanoplot -p
        NanoPlot {params.flags} -o plots/nanoplot -p {wildcards.tag}_alignment_preConsensus_ --bam {input[0]}
        """

rule NanoPlot_alignment:
    input:
        'alignments/{tag}.bam'
    output:
        'plots/nanoplot/{tag, [^\/_]*}_alignment_NanoStats.txt'
    params:
        flags = lambda wildcards: config['nanoplot_flags']
    shell:
        """
        mkdir plots/nanoplot -p
        NanoPlot {params.flags} -o plots/nanoplot -p {wildcards.tag}_alignment_ --bam {input[0]}
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

rule generate_barcode_ref:
    input:
        'alignments/{tag}.bam'
    output:
        'demux/.{tag, [^\/_]*}_generate_barcode_ref.done'
    script:
        'utils/generate_barcode_ref.py'

checkpoint demultiplex:
    input:
        aln = 'alignments/{tag}.bam',
        flag = 'demux/.{tag}_generate_barcode_ref.done'
    output:
        flag = touch('demux/.{tag, [^\/_]*}_demultiplex.done'),
        # checkpoint outputs have the following structure: demux/{tag}_{barcodeGroup}.BAM'
        stats = 'demux/{tag, [^\/_]*}_demux-stats.csv'
    script:
        'utils/demux.py'

rule merge_demux_stats:
    input:
        expand('demux/{tag}_demux-stats.csv', tag=[tag for tag in config['runs'] if config['do_demux'][tag]])
    output:
        'demux-stats.csv'
    run:
        import pandas as pd
        dfs = [pd.read_csv(f) for f in input]
        combined = pd.concat(dfs)
        combined.to_csv(output[0], index=False)

rule index_demuxed:
    input:
        'demux/{tag}_{barcodes}.bam'
    output:
        'demux/{tag, [^\/]*}_{barcodes}.bam.bai'
    shell:
        """
        samtools index {input}
        """

# Mutation analysis will only output AA analysis when a third reference sequence is provided, yielding a dynamic number of output files. Can't use functions in output, so creating a separate
#   rule for which correct input files are only given when AA analysis is not being performed, and giving this rule priority. It's not pretty but it works.
ruleorder: mutation_analysis_NTonly > mutation_analysis

rule mutation_analysis_NTonly:
    input:
        bam = lambda wildcards: 'dummyfilethatshouldneverexist' if config['do_AA_mutation_analysis'][tag] else expand('demux/{tag}_{{barcodes}}.bam', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam',
        bai = lambda wildcards: 'dummyfilethatshouldneverexist' if config['do_AA_mutation_analysis'][tag] else expand('demux/{tag}_{{barcodes}}.bam.bai', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam.bai'
    output:
        expand('mutation_data/{{tag, [^\/_]*}}/{{barcodes, [^\/_]*}}/{{tag}}_{{barcodes}}_{datatype}', datatype = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-muts-frequencies.csv', 'NT-muts-distribution.csv'])
    script:
        'utils/mutation_analysis.py'

rule mutation_analysis:
    input:
        bam = lambda wildcards: expand('demux/{tag}_{{barcodes}}.bam', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam',
        bai = lambda wildcards: expand('demux/{tag}_{{barcodes}}.bam.bai', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam.bai'
    output:
        expand('mutation_data/{{tag, [^\/_]*}}/{{barcodes, [^\/_]*}}/{{tag}}_{{barcodes}}_{datatype}', datatype = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-muts-frequencies.csv', 'NT-muts-distribution.csv', 'AA-muts-frequencies.csv', 'AA-muts-distribution.csv'])
    script:
        'utils/mutation_analysis.py'

def mut_stats_input(wildcards):
    datatypes = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-muts-frequencies.csv', 'NT-muts-distribution.csv']
    if config['do_AA_mutation_analysis'][tag]: datatypes.extend(['AA-muts-frequencies.csv', 'AA-muts-distribution.csv'])
    if config['do_demux'][wildcards.tag]:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split('demultiplex')[0]
        checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
        out = expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{datatype}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, datatype=datatypes)
    else:
        out = expand('mutation_data/{tag}/all/{tag}_all_{datatype}', tag=wildcards.tag, datatype=datatypes)
    assert len(out) > 0, "No demux output files with > minimum count. Pipeline halting."
    return out

rule mut_stats:
	input:
		mut_stats_input
	output:
		'mutation_data/{tag, [^\/]*}/{tag}_mutation-stats.csv'
	script:
		'utils/mutation_statistics.py'

rule merge_mut_stats:
    input:
        expand('mutation_data/{tag}/{tag}_mutation-stats.csv', tag=[tag for tag in config['runs'] if config['do_NT_mutation_analysis'][tag]])
    output:
        'mutation-stats.csv'
    run:
        import pandas as pd
        dfs = [pd.read_csv(f) for f in input]
        combined = pd.concat(dfs)
        combined.to_csv(output[0], index=False)

def dms_view_input(wildcards):
    out = []
    for tag in config['runs']:
        if config['do_demux'][tag]:
            checkpoint_demux_output = checkpoints.demultiplex.get(tag=tag).output[0]
            checkpoint_demux_prefix = checkpoint_demux_output.split('demultiplex')[0]
            checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
            tagFiles = expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_AA-muts-frequencies.csv', tag=tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs)
        else:
            tagFiles = expand('mutation_data/{tag}/{tag}_all_AA-muts-frequencies.csv', tag=tag)
        out.extend(tagFiles)
    assert len(out) > 0, "No demux output files with > minimum count. Pipeline halting."
    return out

# combines all AA mutation frequency data into a table that is compatible with dms-view (https://dms-view.github.io/) which visualizes mutations onto a crystal structure
rule dms_view:
    input:
        dms_view_input
    output:
        'dms-view-table.csv'
    script:
        'utils/frequencies_to_dmsview.py'

rule plot_mutation_spectrum:
    input:
        'mutation_data/{tag}/{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/]*}_mutation-spectra.html'
    script:
        'utils/plot_mutation_spectrum.py'

rule plot_mutation_rate:
    input:
        mutStats = 'mutation-stats.csv',
        timepoints = lambda wildcards: config['timepoints'][wildcards.tag]
    output:
        rate = 'plots/{tag, [^\/]*}_mutation-rates.html',
        rateCSV = 'mutation_data/{tag, [^\/]*}/{tag}_mutation-rates.csv',
        spectrum = 'plots/{tag, [^\/]*}_mutation-rate-spectrum.html',
        spectrumCSV = 'mutation_data/{tag, [^\/]*}/{tag}_mutation-rate-spectrum.csv'
    script:
        'utils/plot_mutation_rate.py'

def plot_mutations_frequencies_input(wildcards):
    if config['do_demux'][wildcards.tag]:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex')[0]
        checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
        return expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{AAorNT}-muts-frequencies.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, AAorNT = wildcards.AAorNT)
    else:
        return expand('mutation_data/{tag}/all/{tag}_all_{AAorNT}-muts-frequencies.csv', tag=wildcards.tag, AAorNT=wildcards.AAorNT)

rule plot_mutations_frequencies:
    input:
        frequencies = plot_mutations_frequencies_input,
        mutStats = 'mutation_data/{tag}/{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/_]*}_{AAorNT, [^\/_]*}-mutations-frequencies.html'
    script:
        'utils/plot_mutations_frequencies.py'

rule plot_mutations_frequencies_barcodeGroup:
    input:
        frequencies = 'mutation_data/{tag}_{barcodes}_{AAorNT}-muts-frequencies.csv',
        mutStats = '{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/_]*}_{barcodes, [^\/_]*}_{AAorNT, [^\/_]*}-muts-frequencies.html'
    script:
        'utils/plot_mutations_frequencies.py'

def plot_mutations_distribution_input(wildcards):
    if config['do_demux'][wildcards.tag]:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex')[0]
        checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
        return expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{AAorNT}-muts-distribution.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, AAorNT = wildcards.AAorNT)
    else:
        return expand('mutation_data/{tag}/all/{tag}_all_{AAorNT}-muts-distribution.csv', tag=wildcards.tag, AAorNT=wildcards.AAorNT)

rule plot_mutations_distribution:
    input:
        dist = plot_mutations_distribution_input
    output:
        'plots/{tag, [^\/_]*}_{AAorNT, [^\/_]*}-mutation-distributions.html'
    script:
        'utils/plot_mutation_distribution.py'

rule plot_mutations_distribution_barcodeGroup:
    input:
        dist = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{AAorNT}-muts-distribution.csv'
    output:
        'plots/{tag, [^\/_]*}_{barcodes, [^\/_]*}_{AAorNT, [^\/_]*}-mutation-distributions.html'
    script:
        'utils/plot_mutation_distribution.py'

rule plot_mutation_diversity:
    input:
        'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv'
    output:
        HamDistPlot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_hamming-distance-distribution.html',
        HamDistCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_hamming-distance-distribution.csv',
        GraphPlot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_diversity-graph.html',
        GraphFile = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_diversity-graph.gexf'
    script:
        'utils/plot_mutation_diversity.py'

def all_diversity_plots_input(wildcards):
    out = []
    checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
    checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex')[0]
    checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
    out.extend( expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{dataType}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, dataType = ['diversity-graph.gexf', 'hamming-distance-distribution.csv']) )
    out.extend( expand('plots/{tag}/{barcodes}/{tag}_{barcodes}_{plotType}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, plotType = ['diversity-graph.html', 'hamming-distance-distribution.html']) )
    return out

rule plot_mutation_diversity_all:
    input:
        all_diversity_plots_input
    output:
        touch('plots/.{tag}_allDiversityPlots.done')

rule plot_pipeline_throughput:
    input:
        initial = lambda wildcards: f'sequences/{wildcards.tag}.fastq.gz',
        UMI_preconsensus_alignment = lambda wildcards: expand('sequences/UMI/{tag}_noConsensus.bam', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_preconsensus_log = lambda wildcards: expand('sequences/UMI/{tag}_noConsensus.log', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_extract = lambda wildcards: expand('sequences/UMI/{tag}_UMI-extract.csv', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_group = lambda wildcards: expand('sequences/UMI/{tag}_UMIgroup-distribution.csv', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_consensus = lambda wildcards: expand('sequences/UMI/{tag}_UMIconsensus.fasta.gz', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        alignment = lambda wildcards: f'alignments/{wildcards.tag}.bam',
        alignment_log = lambda wildcards: f'alignments/{wildcards.tag}.log',
        demux = lambda wildcards: f'demux/{wildcards.tag}_demux-stats.csv' if config['do_demux'][wildcards.tag] else f'sequences/{wildcards.tag}.fastq.gz'
    output:
        plot = 'plots/{tag, [^\/_]*}_pipeline-throughput.html',
        csv = 'maple/{tag, [^\/_]*}_pipeline-throughput.csv'
    script:
        'utils/plot_pipeline_throughput.py'
