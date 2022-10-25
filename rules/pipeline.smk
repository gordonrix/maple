# imports
import os, sys, glob

def get_batch_sequences(runpath, fastq_dir):
    batches_fastqgz, = glob_wildcards("{runpath}/{fastq_dir}/{{id}}.fastq.gz".format(runpath=runpath, fastq_dir=fastq_dir))
    return batches_fastqgz

# get batches of sequences for each provided runname to be merged together
def get_batches(wildcards):
    output = []
    for runname in config['runs'][wildcards.tag]['runname']:

        # search through all directories two directories deep within the minknowDir for the provided runname, and stop once a matching directory is found
        experimentDirs = [d for d in os.listdir(config['minknowDir']) if os.path.isdir(os.path.join(config['minknowDir'],d))]
        expt_sample = False
        for expt in experimentDirs:
            if expt_sample: break # stop once a matching directory is found
            sampleDirs = [d for d in os.listdir(os.path.join(config['minknowDir'], expt)) if os.path.isdir(os.path.join(config['minknowDir'], expt,d))]
            for sample in sampleDirs:
                if runname in os.listdir(os.path.join(config['minknowDir'], expt, sample)):
                    expt_sample = (expt, sample)
                    break

        if expt_sample: # add all batches of sequences to a list to be merged together
            outputs = []
            for fastq_dir in config['fastq_dir'].replace(' ','').split(','):
                batch = get_batch_sequences(os.path.join(config['minknowDir'].rstrip('/'), expt_sample[0], expt_sample[1], runname), fastq_dir)
                outputs.extend(expand("{minknowDir}/{expt}/{sample}/{runname}/{fastq_dir}/{batch}.fastq.gz",
                                        minknowDir = config['minknowDir'].rstrip('/'),
                                        expt = expt_sample[0],
                                        sample = expt_sample[1],
                                        runname = runname,
                                        fastq_dir = fastq_dir,
                                        batch = batch))
        else:
            print('[WARNING] No folders matching the provided runname was found. This is fine if you have already combined the sequences you want to combine but if not then it is not fine and you should refer to the documentation.')
        output.extend(outputs)
    return output

rule merge_paired_end:
    input:
        fwd = lambda wildcards: os.path.join('sequences', 'paired', config['runs'][wildcards.tag]['fwdReads']),
        rvs = lambda wildcards: os.path.join('sequences', 'paired', config['runs'][wildcards.tag]['rvsReads'])
    output:
        merged = temp("sequences/paired/{tag, [^\/_]*}.fastq.gz"),
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

rule basecaller_combine_tag: # combine batches of basecalled reads into a single file
    input:
        lambda wildcards: get_batches(wildcards)
    output:
        temp("sequences/{tag, [^\/_]*}_combined.fastq.gz")
    run:
        with open(output[0], 'wb') as fp_out:
            if len(input)==0:
                raise RuntimeError(f"Basecalled sequence batches not found for tag `{wildcards.tag}`.")
            for f in input:
                with open(f, 'rb') as fp_in:
                    fp_out.write(fp_in.read())

rule move_seqs: # allows for merging batches of sequences or merging paired end reads depending on the tag definition using the above rules
    input:
        lambda wildcards: f'sequences/{wildcards.tag}_combined.fastq.gz' if config['merge_paired_end'][tag]==False else f'sequences/paired/{wildcards.tag}.fastq.gz'
    output:
        'sequences/{tag, [^\/_]*}.fastq.gz'
    shell:
        """
        mv {input} {output}
        """

# default behavior is to use R2C2 as is. The commented out code can use R2C2 just for read splitting, followed by medaka for consensus generation
#  tests showed that medaka produced slightly worse results in most cases. keeping this available for further testing until I am emotionally detached
#  from the time that I spent on this that was likely in vain

rule RCA_consensus:
    input:
        sequence = 'sequences/{tag}.fastq.gz',
        splintRef = lambda wildcards: os.path.join(config['references_directory'], f".{wildcards.tag}_splint.fasta")
    output:
        consensus = 'sequences/RCA/{tag, [^\/_]*}_RCAconsensuses-nofilter.fasta.gz',
        tempOutDir = temp(directory('sequences/tempOutDir-{tag, [^\/_]*}_RCA-consensus')),
        log = 'log/RCA/{tag, [^\/_]*}_RCA.log'
    params:
        peakFinderSettings = lambda wildcards: config['peak_finder_settings'],
        batchSize = lambda wildcards: config['RCA_batch_size']
    threads: workflow.cores
    resources:
        threads = lambda wildcards, threads: threads,
    shell:
        """
        rm -r -f sequences/tempOutDir-{wildcards.tag}_RCA-consensus
        python3 -m C3POa -r {input.sequence} -o sequences/tempOutDir-{wildcards.tag}_RCA-consensus -s {input.splintRef} -n {threads} -co -z --peakFinderSettings {params.peakFinderSettings} --groupSize {params.batchSize}
        mv sequences/tempOutDir-{wildcards.tag}_RCA-consensus/splint/R2C2_Consensus.fasta.gz {output.consensus}
        mkdir log/RCA -p
        mv sequences/tempOutDir-{wildcards.tag}_RCA-consensus/c3poa.log {output.log}
        """

rule filter_RCA_consensus:
    input:
        'sequences/RCA/{tag}_RCAconsensuses-nofilter.fasta.gz'
    output:
        filtered = temp('sequences/RCA/{tag}_RCAconsensuses.fasta.gz'),
        csv = 'log/RCA/{tag}_RCA-log.csv'
    params:
        minimum = lambda wildcards: config['RCA_consensus_minimum'],
        maximum = lambda wildcards: config['RCA_consensus_maximum'],
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    run:
        import gzip
        from Bio import SeqIO
        import os
        import pandas as pd

        # calculate min/max original read length and min/max repeats to use for filtering. repeats reported by C3POa do not include reads at the end,
        #   which can be frequent depending on the method of library preparation
        with open(params.alnRef, 'r') as aln:
            alnSeq = aln.readlines()[1]
        onemerLength = len(alnSeq)
        minLen = len(alnSeq) * params.minimum * 0.9
        maxLen = len(alnSeq) * params.maximum * 1.1
        minRepeats = params.minimum - 2
        maxRepeats = params.maximum - 2

        seqMemoryLimit = 100000 # number of sequences to hold in memory before writing
        seqsOutList = []
        count = 0
        input = input[0]
        for f in [output.filtered, output.filtered[:-3]]:
            if os.path.isfile(f):
                os.remove(f)
        os.makedirs(os.path.split(output.csv)[0], exist_ok=True)
        rows = []
        
        with open(output.filtered[:-3], 'a') as outfile:
            with gzip.open(input, 'rt') as infile:
                for record in SeqIO.parse(infile, 'fasta'):
                    readName, avgQual, originalLen, numRepeats, subreadLen = record.id.split('_')
                    passedConsensus = 0
                    if ( minRepeats <= int(numRepeats) <= maxRepeats ) and ( minLen <= float(originalLen) <= maxLen):
                        passedConsensus = 1
                        seqsOutList.append(record)
                        count += 1
                    rows.append([readName, avgQual, originalLen, numRepeats, subreadLen, passedConsensus])
                    if count > seqMemoryLimit:
                        SeqIO.write(seqsOutList, outfile, 'fasta-2line')
                        seqsOutList = []
                        count = 0
            SeqIO.write(seqsOutList, outfile, 'fasta-2line')
        os.system(f'gzip {output.filtered[:-3]}')

        pd.DataFrame(rows, columns=['read_ID', 'average_quality_score', 'original_read_length', 'number_of_complete_repeats', 'subread length', 'pass(1)/fail(0)']).to_csv(output.csv, index=False)

# rule RCA_split_reads:
#     input:
#         sequence = 'sequences/{tag}.fastq.gz',
#         splintRef = lambda wildcards: os.path.join(config['references_directory'], f".{wildcards.tag}_splint.fasta")
#     output:
#         subreads = temp('sequences/RCA/{tag, [^\/_]*}_RCA-subreads.fastq.gz'),
#         log = 'sequences/RCA/{tag, [^\/_]*}_c3poa.log'
#     params:
#         peakFinderSettings = lambda wildcards: config['peak_finder_settings'],
#         minimum = lambda wildcards: config['RCA_consensus_minimum'],
#         maximum = lambda wildcards: config['RCA_consensus_maximum']
#     threads: workflow.cores
#     resources:
#         threads = lambda wildcards, threads: threads
#     shell:
#         """
#         rm -r -f sequences/tempOutDir-{wildcards.tag}_RCA
#         python3 -m C3POa -r {input.sequence} -o sequences/tempOutDir-{wildcards.tag}_RCA -s {input.splintRef} -n {threads} -co -z --peakFinderSettings {params.peakFinderSettings} --groupSize 1000 --minimum {params.minimum} --maximum {params.maximum}
#         mv sequences/tempOutDir-{wildcards.tag}_RCA/splint/R2C2_Subreads.fastq.gz {output.subreads}
#         mv sequences/tempOutDir-{wildcards.tag}_RCA/c3poa.log {output.log}
#         rm -r sequences/tempOutDir-{wildcards.tag}_RCA
#         """

# RCAbatchesList = [str(x) for x in range(0,config['RCA_medaka_batches'])]
# rule RCA_split_fasta:
#     input:
#         subreads = 'sequences/RCA/{tag}_RCA-subreads.fastq.gz'
#     output:
#         fastqs = temp(expand('sequences/RCA/{{tag, [^\/_]*}}-temp/batch{x}.fasta', x=RCAbatchesList))
#     params:
#         batches = config['RCA_medaka_batches']
#     script:
#         'utils/RCA_split_fasta.py'

# rule RCA_consensus:
#     input:
#         fastq = 'sequences/RCA/{tag}-temp/{batch}.fasta',
#         alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
#     output:
#         outDir = temp(directory('sequences/RCA/{tag, [^\/_]*}-temp/{batch, [^\/_]*}')),
#         consensus = temp('sequences/RCA/{tag, [^\/_]*}-temp/{batch, [^\/_]*}/consensus.fasta')
#     threads: lambda wildcards: config['threads_medaka']
#     resources:
#         threads = lambda wildcards, threads: threads,
#         mem_mb = 1000
#     params:
#         model = lambda wildcards: config['medaka_model'],
#         flags = lambda wildcards: config['medaka_flags'],
#         maxmem = lambda wildcards, threads, resources: resources.mem_mb * 1024
#     shell:
#         """
#         rm -rf {output.outDir}
#         medaka maple_smolecule --threads {threads} --model {params.model} {params.flags} {output.outDir} {input.alnRef} {input.fastq}
#         """

# rule RCA_merge_consensus_seqs:
#     input:
#         expand('sequences/RCA/{{tag}}-temp/batch{x}/consensus.fasta', x = RCAbatchesList)
#     output:
#         seqs = 'sequences/RCA/{tag, [^\/_]*}_RCAconsensuses.fasta.gz',
#         log = 'sequences/RCA/{tag, [^\/_]*}_RCAconsensuses.log'
#     run:
#         with open(output.seqs[:-3], 'w') as fp_out, open(output.log, 'w') as log_out:
#             for f in input:
#                 with open(f, 'r') as fp_in:
#                     fp_out.write(fp_in.read())
#         os.system(f'gzip {output.seqs[:-3]}')

rule UMI_minimap2:
    input:
        sequence = lambda wildcards: 'sequences/RCA/{tag, [^\/_]*}_RCAconsensuses.fasta.gz' if config['do_RCA_consensus'][wildcards.tag] else 'sequences/{tag}.fastq.gz',
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        aln = pipe("sequences/UMI/{tag, [^\/_]*}_noConsensus.sam"),
        log = "sequences/UMI/{tag, [^\/_]*}_noConsensus.log"
    threads: config['threads_alignment']
    group: "minimap2"
    params:
        flags = lambda wildcards: config['alignment_minimap2_flags'] if type(config['alignment_minimap2_flags'])==str else config['alignment_minimap2_flags'][wildcards.tag]
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    shell:
        """
        {config[bin_singularity][minimap2]} -t {threads} {params.flags} {input.alnRef} {input.sequence} 1>> {output.aln} 2> >(tee {output.log} >&2)
        if [ $(grep 'ERROR' {output.log} | wc -l) -gt 0 ]; then exit 1; fi
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
        bam = temp('sequences/UMI/{tag, [^\/_]*}_UMIgroup.bam'),
        log = temp('sequences/UMI/{tag, [^\/_]*}_UMIgroup-log.tsv')
    params:
        UMI_mismatches = lambda wildcards: config['UMI_mismatches']
    shell:
        """
        umi_tools group -I {input.bam} --group-out={output.log} --output-bam --per-gene --gene-tag=GN --edit-distance-threshold {params.UMI_mismatches} -S {output.bam}
        """

rule plot_UMI_group:
    input:
        'sequences/UMI/{tag}_UMIgroup-log.tsv'
    output:
        csv = 'sequences/UMI/{tag, [^\/_]*}_UMIgroup-distribution.csv',
        plot = 'plots/{tag, [^\/_]*}_UMIgroup-distribution.html'
    script:
        'utils/plot_UMI_groups_distribution.py'

UMIbatchesList = [str(x) for x in range(0,config['UMI_medaka_batches'])]
rule split_BAMs_to_fasta:
    input:
        grouped = 'sequences/UMI/{tag}_UMIgroup.bam',
        log = 'sequences/UMI/{tag}_UMIgroup-log.tsv'
    output:
        fastas = expand('sequences/UMI/{{tag, [^\/_]*}}-temp/batch{x}.fasta', x=UMIbatchesList)
    params:
        batches = lambda wildcards: config['UMI_medaka_batches'],
        minimum = lambda wildcards: config['UMI_consensus_minimum'],
        maximum = lambda wildcards: config['UMI_consensus_maximum']
    script:
        'utils/UMI_splitBAMs.py'

# use medaka to generate consensus sequences if either min or max reads/UMI not set to 1, otherwise can just merge the sequences as they have been deduplicated by the split_BAMs_to_fasta rule
if config['UMI_consensus_minimum'] == config['UMI_consensus_maximum'] == 1:

    rule UMI_merge_deduplicated_seqs:
        input:
            expand('sequences/UMI/{{tag}}-temp/batch{x}.fasta', x = UMIbatchesList)
        output:
            seqs = 'sequences/UMI/{tag, [^\/_]*}_UMIconsensuses.fasta.gz',
            log = 'sequences/UMI/{tag, [^\/_]*}_UMIconsensuses.log'
        run:
            with open(output.seqs[:-3], 'w') as fp_out, open(output.log, 'w') as log_out:
                for f in input:
                    with open(f, 'r') as fp_in:
                        fp_out.write(fp_in.read())
            os.system(f'gzip {output.seqs[:-3]}')

else:

    rule UMI_consensus:
        input:
            fasta = 'sequences/UMI/{tag}-temp/{batch}.fasta',
            alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
        output:
            outDir = temp(directory('sequences/UMI/{tag, [^\/_]*}-temp/{batch, [^\/_]*}')),
            consensus = temp('sequences/UMI/{tag, [^\/_]*}-temp/{batch, [^\/_]*}_consensus.fasta')
        params:
            depth = lambda wildcards: config['UMI_consensus_minimum'],
            model = lambda wildcards: config['medaka_model'],
            flags = lambda wildcards: config['medaka_flags']
        threads: workflow.cores
        resources:
            threads = lambda wildcards, threads: threads
        shell:
            """
            rm -rf {output.outDir}
            medaka maple_smolecule --threads {threads} --model {params.model} {params.flags} --depth {params.depth} {output.outDir} {input.alnRef} {input.fasta}
            mv {output.outDir}/consensus.fasta {output.consensus}
            """

    rule UMI_merge_consensus_seqs:
        input:
            expand('sequences/UMI/{{tag}}-temp/batch{x}_consensus.fasta', x = UMIbatchesList)
        output:
            seqs = 'sequences/UMI/{tag, [^\/_]*}_UMIconsensuses.fasta.gz'
        run:
            with open(output.seqs[:-3], 'w') as fp_out:
                for f in input:
                    with open(f, 'r') as fp_in:
                        fp_out.write(fp_in.read())
            os.system(f'gzip {output.seqs[:-3]}')

def alignment_sequence_input(wildcards):
    if config['do_UMI_analysis'][wildcards.tag]:
        return expand('sequences/UMI/{tag}_UMIconsensuses.fasta.gz', tag=config['consensusCopyDict'][wildcards.tag])
    elif config['do_RCA_consensus'][wildcards.tag]:
        return 'sequences/RCA/{tag}_RCAconsensuses.fasta.gz'
    else:
        return 'sequences/{tag}.fastq.gz'

rule minimap2:
    input:
        sequence = alignment_sequence_input,
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        aln = pipe("alignments/{tag, [^\/_]*}.sam"),
        log = "alignments/{tag, [^\/_]*}.log"
    params:
        flags = lambda wildcards: config['alignment_minimap2_flags'] if type(config['alignment_minimap2_flags'])==str else config['alignment_minimap2_flags'][wildcards.tag]
    threads: config['threads_alignment']
    group: "minimap2"
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt * config['runtime']['minimap2'])   # 60 min / 16 threads
    shell:
        """
        {config[bin_singularity][minimap2]} -t {threads} {params.flags} {input.alnRef} {input.sequence} 1>> {output.aln} 2> >(tee {output.log} >&2)
        if [ $(grep 'ERROR' {output.log} | wc -l) -gt 0 ]; then exit 1; fi
        """

# sam to bam conversion
rule aligner_sam2bam:
    input:
        sam = "alignments/{tag}.sam"
    output:
        bam = temp("alignments/{tag, [^\/_]*}.bam"),
        bai = temp("alignments/{tag, [^\/_]*}.bam.bai")
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
        bai = 'alignments/{tag}.bam.bai',
        flag = 'demux/.{tag}_generate_barcode_ref.done'
    output:
        flag = touch('demux/.{tag, [^\/_]*}_demultiplex.done'),
        # checkpoint outputs have the following structure: demux/{tag}_{barcodeGroup}.bam'
        stats = 'demux/{tag, [^\/_]*}_demux-stats.csv'
    params:
        barcodeInfo = lambda wildcards: config['runs'][wildcards.tag]['barcodeInfo'],
        barcodeGroups = lambda wildcards: config['runs'][wildcards.tag].get('barcodeGroups', False)
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

def ma_NTonly_input(wildcards):
    if config['do_AA_mutation_analysis'][wildcards.tag]:
        return {'bam':'dummyfilethatshouldneverexist','bai':'dummyfilethatshouldneverexist'}
    else:
        if config['do_demux'][wildcards.tag]:
            return {'bam':expand('demux/{tag}_{{barcodes}}.bam', tag=wildcards.tag), 'bai':expand('demux/{tag}_{{barcodes}}.bam.bai', tag=wildcards.tag)}
        else:
            return {'bam':f'alignments/{wildcards.tag}.bam', 'bai':f'alignments/{wildcards.tag}.bam.bai'}

rule mutation_analysis_NTonly:
    input:
        unpack(ma_NTonly_input)
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
    if config['do_AA_mutation_analysis'][wildcards.tag]: datatypes.extend(['AA-muts-frequencies.csv', 'AA-muts-distribution.csv'])
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

ruleorder: mutation_diversity_NTonly > mutation_diversity

rule mutation_diversity_NTonly:
    input:
        lambda wildcards: 'dummyfilethatshouldneverexist' if config['do_AA_mutation_analysis'][wildcards.tag] else 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv'
    output:
        ntHamDistCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_NT-hamming-distance-distribution.csv',
        edges = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_edges.csv'
    params:
        downsample = lambda wildcards: config.get('diversity_plot_downsample', False)
    script:
        'utils/mutation_diversity.py'

rule mutation_diversity:
    input:
        'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv'
    output:
        ntHamDistCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_NT-hamming-distance-distribution.csv',
        aaHamDistCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_AA-hamming-distance-distribution.csv',
        edges = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_edges.csv'
    params:
        downsample = lambda wildcards: config.get('diversity_plot_downsample', False)
    script:
        'utils/mutation_diversity.py'

rule plot_mutation_diversity:
    input:
        genotypes = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv',
        edges = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_edges.csv',
        ntHamDistCSV = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_NT-hamming-distance-distribution.csv'
    output:
        ntHamDistPlot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_NT-hamming-distance-distribution.html',
        GraphPlot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_diversity-graph.html',
        GraphFile = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_diversity-graph.gexf'
    params:
        edgeLimit = lambda wildcards: config.get('diversity_plot_hamming_distance_edge_limit', False),
        xMax = lambda wildcards: config['hamming_distance_distribution_plot_x_max'],
        nodeSize = lambda wildcards: config['force_directed_plot_node_size'],
        nodeColor = lambda wildcards: config['force_directed_plot_node_color'],
        downsample = lambda wildcards: config.get('diversity_plot_downsample', False)
    script:
        'utils/plot_mutation_diversity.py'

rule plot_AA_mutation_diversity:
    input:
        aaHamDistCSV = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_AA-hamming-distance-distribution.csv'
    output:
        aaHamDistPlot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_AA-hamming-distance-distribution.html'
    params:
        xMax = lambda wildcards: config['hamming_distance_distribution_plot_x_max'],
        downsample = lambda wildcards: config.get('diversity_plot_downsample', False)
    script:
        'utils/plot_mutation_diversity.py'

def all_diversity_plots_input(wildcards):
    out = []
    checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
    checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex')[0]
    checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
    out.extend( expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{dataType}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, dataType = ['diversity-graph.gexf', 'NT-hamming-distance-distribution.csv']) )
    out.extend( expand('plots/{tag}/{barcodes}/{tag}_{barcodes}_{plotType}', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, plotType = ['diversity-graph.html', 'NT-hamming-distance-distribution.html']) )
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
        UMI_preconsensus_index = lambda wildcards: expand('sequences/UMI/{tag}_noConsensus.bam.bai', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_preconsensus_log = lambda wildcards: expand('sequences/UMI/{tag}_noConsensus.log', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_extract = lambda wildcards: expand('sequences/UMI/{tag}_UMI-extract.csv', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_group = lambda wildcards: expand('sequences/UMI/{tag}_UMIgroup-distribution.csv', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        UMI_consensus = lambda wildcards: expand('sequences/UMI/{tag}_UMIconsensuses.fasta.gz', tag=config['consensusCopyDict'][wildcards.tag])[0] if config['do_UMI_analysis'][wildcards.tag]==True else f'sequences/{wildcards.tag}.fastq.gz',
        alignment = lambda wildcards: f'alignments/{wildcards.tag}.bam',
        alignment_index = lambda wildcards: f'alignments/{wildcards.tag}.bam.bai',
        alignment_log = lambda wildcards: f'alignments/{wildcards.tag}.log',
        demux = lambda wildcards: f'demux/{wildcards.tag}_demux-stats.csv' if config['do_demux'][wildcards.tag] else f'sequences/{wildcards.tag}.fastq.gz'
    output:
        plot = 'plots/{tag, [^\/_]*}_pipeline-throughput.html',
        csv = 'maple/{tag, [^\/_]*}_pipeline-throughput.csv'
    script:
        'utils/plot_pipeline_throughput.py'
