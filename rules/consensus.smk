#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Rules for generating more accurate consensus sequences from sequencing data that contains multiple
#                   reads of the same original DNA sequence
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

# Rolling circle amplification (RCA)-based consensus using C3POa
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
        log = 'log/RCA/{tag}_RCA-log.csv'
    params:
        minimum = lambda wildcards: config['RCA_consensus_minimum'],
        maximum = lambda wildcards: config['RCA_consensus_maximum'],
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    run:
        import gzip
        from Bio import SeqIO
        import os
        import pandas as pd

        # calculate min/max original read length and min/max repeats to use for filtering. repeats reported by C3POa do not include
        #   subreads at the ends of the read that are almost complete, which can be common and valuable depending on the method of library preparation
        with open(params.alnRef, 'r') as aln:
            alnSeq = aln.readlines()[1]
        onemerLength = len(alnSeq)
        minLen = len(alnSeq) * 0.8
        maxLen = len(alnSeq) * 1.2
        minRepeats = params.minimum - 2
        maxRepeats = params.maximum - 2

        seqMemoryLimit = 100000 # number of sequences to hold in memory before writing
        seqsOutList = []
        count = 0
        input = input[0]
        for f in [output.filtered, output.filtered[:-3]]:
            if os.path.isfile(f):
                os.remove(f)
        os.makedirs(os.path.split(output.log)[0], exist_ok=True)
        rows = []
        
        with open(output.filtered[:-3], 'a') as outfile:
            with gzip.open(input, 'rt') as infile:
                for record in SeqIO.parse(infile, 'fasta'):
                    readName, avgQual, originalLen, numRepeats, subreadLen = record.id.split('_')
                    passedConsensus = 'fail'
                    if ( minRepeats <= int(numRepeats) <= maxRepeats ) and ( minLen <= float(subreadLen) <= maxLen):
                        passedConsensus = 'pass'
                        seqsOutList.append(record)
                        count += 1
                    rows.append([readName, avgQual, originalLen, numRepeats, subreadLen, int(originalLen)/int(subreadLen), passedConsensus])
                    if count > seqMemoryLimit:
                        SeqIO.write(seqsOutList, outfile, 'fasta-2line')
                        seqsOutList = []
                        count = 0
            SeqIO.write(seqsOutList, outfile, 'fasta-2line')
        os.system(f'gzip {output.filtered[:-3]}')

        pd.DataFrame(rows, columns=['read_ID', 'average_quality_score', 'original_read_length', 'number_of_complete_repeats', 'subread_length', 'read_length ÷ subread_length', 'pass/fail']).to_csv(output.log, index=False)

rule plot_RCA_consensus:
    input:
        log = 'log/RCA/{tag}_RCA-log.csv'
    output:
        plot = 'plots/{tag, [^\/_]*}_RCA-distribution.html'
    run:
        import pandas as pd
        import holoviews as hv
        hv.extension('bokeh')

        data = pd.read_csv(input.log)[['subread_length','read_length ÷ subread_length','pass/fail']]

        plot = hv.Points(data, vdims=['pass/fail']).hist(dimension=['subread_length','read_length ÷ subread_length'], num_bins=50).opts(
            hv.opts.Points(width=750,height=750,fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}, size=5, alpha=0.3, color='pass/fail', cmap={'pass':'#3296FA','fail':'#FA3232'})).opts(
            hv.opts.Histogram(tools=['hover'],fontsize={'title':16,'labels':14,'xticks':10,'yticks':10},color='grey'))

        hv.save(plot, output.plot, backend='bokeh')

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
