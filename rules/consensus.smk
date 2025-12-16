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
        splintRef = lambda wildcards: os.path.join(config['metadata'], f".{wildcards.tag}_splint.fasta")
    output:
        consensus = 'sequences/RCA/{tag, [^\/_]*}_RCAconsensuses-nofilter.fasta.gz',
        tempOutDir = temp(directory('sequences/tempOutDir-{tag, [^\/_]*}_RCA-consensus')),
        log = 'log/RCA/{tag, [^\/_]*}_RCA.log'
    params:
        peakFinderSettings = lambda wildcards: config.get('peak_finder_settings', '23,3,27,2'),
        batchSize = lambda wildcards: config['RCA_batch_size']
    threads: workflow.cores
    resources:
        threads = lambda wildcards, threads: threads,
    shell:
        """
        rm -r -f sequences/tempOutDir-{wildcards.tag}_RCA-consensus
        python3 -m C3POa -r {input.sequence} -o sequences/tempOutDir-{wildcards.tag}_RCA-consensus -s {input.splintRef} -n {threads} -co -z --peakFinderSettings {params.peakFinderSettings}
        mv sequences/tempOutDir-{wildcards.tag}_RCA-consensus/splint/R2C2_Consensus.fasta.gz {output.consensus}
        mkdir -p log/RCA
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
            hv.opts.Points(width=500,height=500,fontsize={'title':16,'labels':14,'xticks':10,'yticks':10}, size=5, alpha=0.3, color='pass/fail', cmap={'pass':'#3296FA','fail':'#FA3232'})).opts(
            hv.opts.Histogram(tools=['hover'],fontsize={'title':16,'labels':14,'xticks':10,'yticks':10},color='grey'))

        hv.save(plot, output.plot, backend='bokeh')

rule UMI_align:
    input:
        sequence = lambda wildcards: 'sequences/RCA/{tag, [^\/_]*}_RCAconsensuses.fasta.gz' if config['do_RCA_consensus'][wildcards.tag] else 'sequences/{tag}.fastq.gz',
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        aln = temp("sequences/UMI/{tag, [^\/_]*}_pre-consensus.sam"),
        log = "sequences/UMI/{tag, [^\/_]*}_pre-consensus.log"
    threads: config['threads_alignment']
    params:
        flags = lambda wildcards: config['alignment_minimap2_flags'] if type(config['alignment_minimap2_flags'])==str else config['alignment_minimap2_flags'][wildcards.tag]
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, threads, attempt: int((1.0 + (0.2 * (attempt - 1))) * (config['memory']['minimap2'][0] + config['memory']['minimap2'][1] * threads)),
        time_min = lambda wildcards, threads, attempt: int((960 / threads) * attempt)   # 60 min / 16 threads
    shell:
        """
        minimap2 -t {threads} {params.flags} {input.alnRef} {input.sequence} 1>> {output.aln} 2> >(tee {output.log} >&2)
        if [ $(grep 'ERROR' {output.log} | wc -l) -gt 0 ]; then exit 1; fi
        """

# sam to bam conversion
rule UMI_aligner_sam2bam:
    input:
        sam = "sequences/UMI/{tag}_pre-consensus.sam"
    output:
        bam = temp("sequences/UMI/{tag, [^\/_]*}_pre-consensus.bam"),
        index = temp("sequences/UMI/{tag, [^\/_]*}_pre-consensus.bam.bai")
    shadow: "minimal"
    threads: 1
    resources:
        threads = lambda wildcards, threads: threads,
        mem_mb = lambda wildcards, attempt: int((1.0 + (0.2 * (attempt - 1))) * 5000)
    # singularity:
    #     config['singularity_images']['alignment']
    shell:
        """
        samtools view -b {input.sam} | samtools sort -m 4G > {output.bam}
        samtools index {output.bam}
        """

# Unified UMI processing: extract UMIs from aligned BAM, group by (reference, UMI), filter, and split into batches
UMIbatchesList = [str(x) for x in range(0, config['UMI_medaka_batches'])]
rule UMI_process:
    input:
        bam = 'sequences/UMI/{tag}_pre-consensus.bam',
        index = "sequences/UMI/{tag}_pre-consensus.bam.bai",
        references = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        bams = temp(expand("sequences/UMI/{{tag, [^\/_]*}}-temp/batch{x}.bam", x=UMIbatchesList)),
        bais = temp(expand("sequences/UMI/{{tag, [^\/_]*}}-temp/batch{x}.bam.bai", x=UMIbatchesList)),
        extract_log = temp('sequences/UMI/{tag, [^\/_]*}_UMI-extract.csv'),
        groups_log = temp("sequences/UMI/{tag, [^\/_]*}_UMI-groups-log.csv")
    params:
        UMI_contexts = lambda wildcards: config['runs'][wildcards.tag]['UMI_contexts'],
        UMI_mismatches = lambda wildcards: config['UMI_mismatches'],
        batches = lambda wildcards: config['UMI_medaka_batches'],
        minimum = lambda wildcards: config['UMI_consensus_minimum'],
        maximum = lambda wildcards: config['UMI_consensus_maximum'],
        output_format = 'bam'
    script:
        "utils/UMI_process.py"

# # OLD: extract UMI sequences from the bam file and save both original sequences and UMI sequences to fastq files
# rule UMI_extract:
#     input:
#         bam = 'sequences/UMI/{tag}_noConsensus.bam',
#         index = "sequences/UMI/{tag}_noConsensus.bam.bai"
#     output:
#         UMIs = temp('sequences/UMI/{tag, [^\/_]*}_UMIs.fastq'),
#         sequences = temp('sequences/UMI/{tag, [^\/_]*}_sequences.fastq'),
#         log = temp('sequences/UMI/{tag, [^\/_]*}_UMI-extract.csv')
#     params:
#         reference = lambda wildcards: config['runs'][wildcards.tag]['reference'],
#         UMI_contexts = lambda wildcards: ",".join(config['runs'][wildcards.tag]['UMI_contexts']),
#         mode = "fastq"
#     script:
#         "utils/UMI_extract.py"


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
        outDF = pd.DataFrame(columns=['tag','reference_name','success','failure']+UMIcols)

        for f in input:
            df = pd.read_csv(f, dtype={'umi':str})
            tag = f.split('/')[-1].split('_')[0]
            df['tag'] = tag

            # Group by tag and reference_name, sum the numeric columns
            sum_cols = ['success', 'failure'] + [col for col in df.columns if col.startswith('umi_') and col.endswith('_failure')]
            dfSum = df.groupby(['tag', 'reference_name'])[sum_cols].sum().reset_index()

            # Add missing UMI columns if needed
            for col in UMIcols:
                if col not in dfSum.columns:
                    dfSum[col] = 'NA'

            # Sort this tag's data by success (descending) then reference_name
            dfSum = dfSum.sort_values(by=['success', 'reference_name'], ascending=[False, True])

            outDF = pd.concat([outDF, dfSum], ignore_index=True)

        outDF.to_csv(output[0], index=False)

# # OLD: Group UMIs using UMICollapse
# rule UMI_group:
#     input:
#         UMIs = 'sequences/UMI/{tag}_UMIs.fastq'
#     output:
#         UMIs_grouped = temp('sequences/UMI/{tag, [^\/_]*}_UMIs-grouped.fastq')
#     params:
#         UMI_mismatches = lambda wildcards: config['UMI_mismatches']
#     shell:
#         """
#         umicollapse fastq -i {input.UMIs} -o {output.UMIs_grouped} --tag -k {params.UMI_mismatches}
#         """

# # OLD: transfer the cluster ID from the UMI file to the corresponding sequence
# rule transfer_cluster_id:
#     input:
#         umi_fastq = 'sequences/UMI/{tag}_UMIs-grouped.fastq',
#         seq_fastq = 'sequences/UMI/{tag}_sequences.fastq'
#     output:
#         fastq = temp('sequences/UMI/{tag, [^\/_]*}_sequences-grouped.fastq')
#     script:
#         "utils/UMI_transfer_id.py"

# # OLD: Generate UMI groups log
# rule UMI_groups_log:
#     input:
#         fastq = "sequences/UMI/{tag}_UMIs-grouped.fastq"
#     output:
#         csv = temp("sequences/UMI/{tag, [^\/_]*}_UMI-groups-log.csv")
#     script:
#         "utils/UMI_groups_log.py"

rule UMI_group_distribution:
    input:
        csv = 'sequences/UMI/{tag}_UMI-groups-log.csv'
    output:
        csv = 'sequences/UMI/{tag, [^\/_]*}_all_UMIgroup-distribution.csv'
    run:
        import pandas as pd
        import numpy as np
        from utils.common import dist_to_DF

        df = pd.read_csv(input.csv)

        # Process each reference separately
        all_dfs = []
        for ref_name in df['reference_name'].unique():
            ref_df = df[df['reference_name'] == ref_name]
            df_subset = ref_df[['final_umi_count', 'unique_id']].dropna()
            df_unique = df_subset.drop_duplicates()
            UMIcounts = df_unique['final_umi_count'].to_numpy().astype(int)
            dist = np.bincount(UMIcounts)
            outDF = dist_to_DF(dist, 'subreads', 'UMI groups')
            outDF = outDF.loc[1:, :] # skip first row because there are never 0 subreads
            outDF['reference_name'] = ref_name
            all_dfs.append(outDF)

        # Concatenate all reference distributions
        final_df = pd.concat(all_dfs, ignore_index=True)
        final_df.to_csv(output.csv, index=False)

rule plot_distribution_UMI_group:
    input:
        expand('sequences/UMI/{tag}_all_UMIgroup-distribution.csv', tag = list(set( [config['consensusCopyDict'][str(t)] for t in config['runs'] if config['do_UMI_analysis'][t]] )))
    output:
        plot = 'plots/UMI-group-distribution.html'
    params:
        labels = list(set( [config['consensusCopyDict'][str(t)] for t in config['runs'] if config['do_UMI_analysis'][t]] )),
        title = 'all tags',
        legend_label = 'tag',
        background = False,
        raw = True,
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        x_max = lambda wildcards: config.get('distribution_x_range', False),
        y_max = lambda wildcards: config.get('distribution_y_range', False),
        log_x = lambda wildcards: config.get('UMI_distribution_log_x', False),
        log_y = lambda wildcards: config.get('UMI_distribution_log_y', False),
        split_by_reference = lambda wildcards: config.get('split_by_reference', False)
    script:
        'utils/plot_distribution.py'

# # OLD: Split sequences into batches for consensus
# rule split_FASTQs_to_fasta:
#     input:
#         fastq="sequences/UMI/{tag}_sequences-grouped.fastq",
#         log="sequences/UMI/{tag}_UMI-groups-log.csv"
#     output:
#         fastas = temp(expand("sequences/UMI/{{tag, [^\/_]*}}-temp/batch{x}.fasta", x=UMIbatchesList))
#     params:
#         batches = lambda wildcards: config['UMI_medaka_batches'],
#         minimum = lambda wildcards: config['UMI_consensus_minimum'],
#         maximum = lambda wildcards: config['UMI_consensus_maximum']
#     script:
#         "utils/UMI_split_fastqs.py"

# use medaka to generate consensus sequences if either min or max reads/UMI not set to 1, otherwise can just merge the sequences as they have been deduplicated by the UMI_process rule
if config['UMI_consensus_minimum'] == config['UMI_consensus_maximum'] == 1:

    rule UMI_merge_deduplicated_seqs:
        input:
            bams = expand('sequences/UMI/{{tag}}-temp/batch{x}.bam', x = UMIbatchesList)
        output:
            seqs = 'sequences/UMI/{tag, [^\/_]*}_UMIconsensuses.fasta.gz'
        run:
            import pysam
            with open(output.seqs[:-3], 'w') as fp_out:
                for bam_file in input.bams:
                    with pysam.AlignmentFile(bam_file, 'rb', check_sq=False) as bam:
                        for read in bam:
                            fp_out.write(f'>{read.query_name}\n{read.query_sequence}\n')
            os.system(f'gzip {output.seqs[:-3]}')

else:

    rule UMI_consensus:
        input:
            bam = 'sequences/UMI/{tag}-temp/batch{batch}.bam',
            bai = 'sequences/UMI/{tag}-temp/batch{batch}.bam.bai',
            references = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
        output:
            outDir = temp(directory('sequences/UMI/{tag, [^\/_]*}-temp/batch{batch,\\d+}')),
            consensus = temp('sequences/UMI/{tag, [^\/_]*}-temp/batch{batch,\\d+}_consensus.fasta')
        params:
            depth = lambda wildcards: config['UMI_consensus_minimum'],
            model = lambda wildcards: config['medaka_model'],
            flags = lambda wildcards: config['medaka_flags']
        threads: workflow.cores if config.get('medaka_use_gpu', True) else 1
        resources:
            threads = lambda wildcards, threads: threads
        shell:
            """
            echo "Running medaka consensus with {threads} threads"
            rm -rf {output.outDir}
            medaka maple_smolecule --threads {threads} --model {params.model} {params.flags} --depth {params.depth} {output.outDir} {input.references} {input.bam}
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
