#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Rules for various demultiplexing, analysis, and plotting steps
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

# imports
import os, sys, glob
import pandas as pd
from utils.common import sort_barcodes, dashboard_input, process_timepoints_csv

# Prefer timepoint-specific rule over generic rule for timepoint files
ruleorder: reduce_timepoint_genotypes_dimensions > reduce_genotypes_dimensions

def get_dashboard_datasets_df():
    """Generate the dashboard datasets dataframe based on config"""
    rows = []

    # If timepoints exist, only add timepoints
    if 'timepointsInfo' in config:
        if config.get('do_enrichment', False):
            # With enrichment: one entry per group (aggregated across replicates)
            seen_groups = set()
            for sample_label, sample_info in config['timepointsInfo'].items():
                group = sample_info['group']
                if group in seen_groups:
                    continue
                seen_groups.add(group)

                enrichment_mode = config.get('enrichment_type', 'genotype')
                genotypes_csv = f'mutation_data/timepoints/{group}_merged-group_genotypes-reduced-dimensions-{enrichment_mode}-enrichment-mean.csv'
                reference_csv = sample_info.get('reference_csv', '')
                exclude_indels = str(not config.get('analyze_seqs_with_indels', True))

                rows.append({
                    'genotypes_csv': genotypes_csv,
                    'reference_csv': reference_csv,
                    'exclude_indels': exclude_indels
                })
        else:
            # Without enrichment: one entry per sample
            for sample_label, sample_info in config['timepointsInfo'].items():
                genotypes_csv = f'mutation_data/timepoints/{sample_label}_merged-timepoint_genotypes-reduced-dimensions.csv'
                reference_csv = sample_info.get('reference_csv', '')
                exclude_indels = str(not config.get('analyze_seqs_with_indels', True))

                rows.append({
                    'genotypes_csv': genotypes_csv,
                    'reference_csv': reference_csv,
                    'exclude_indels': exclude_indels
                })
    # Otherwise add individual tags
    else:
        for tag in config['runs']:
            genotypes_csv = f'mutation_data/{tag}/{tag}_genotypes-reduced-dimensions.csv'
            if 'enrichment' in config['runs'][tag]:
                genotypes_csv = genotypes_csv[:-4] + '-enrichment.csv'

            reference_csv = config['runs'][tag].get('reference_csv', '')
            exclude_indels = str(not config.get('analyze_seqs_with_indels', True))

            rows.append({
                'genotypes_csv': genotypes_csv,
                'reference_csv': reference_csv,
                'exclude_indels': exclude_indels
            })

    return pd.DataFrame(rows)

def alignment_sequence_input(wildcards):
    if config['do_UMI_analysis'][wildcards.tag]:
        return expand('sequences/UMI/{tag}_UMIconsensuses.fasta.gz', tag=config['consensusCopyDict'][wildcards.tag])
    elif config['do_RCA_consensus'][wildcards.tag]:
        return 'sequences/RCA/{tag}_RCAconsensuses.fasta.gz'
    else:
        return 'sequences/{tag}.fastq.gz'

rule align:
    input:
        sequence = alignment_sequence_input,
        alnRef = lambda wildcards: config['runs'][wildcards.tag]['reference_aln']
    output:
        aln = temp("alignments/{tag, [^\/_]*}.sam"),
        log = "alignments/{tag, [^\/_]*}.log"
    params:
        flags = lambda wildcards: config['alignment_minimap2_flags'] if type(config['alignment_minimap2_flags'])==str else config['alignment_minimap2_flags'][wildcards.tag]
    threads: config['threads_alignment'] if config['threads_alignment']<(workflow.cores-1) else max(workflow.cores-1,1)
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
    shell:
        """
        samtools view -b {input.sam} | samtools sort -m 4G > {output.bam}
        samtools index {output.bam}
        """

# Use NanoPlot package to generate and organize plots that describe sequencing data and alignments
# use a try except instead of requiring the nanopore option in config
rule NanoPlot_fastq:
    input:
        'sequences/{tag}.fastq.gz'
    output:
        'plots/nanoplot/{tag, [^\/_]*}_fastq_NanoStats.txt'
    params:
        flags = lambda wildcards: config['nanoplot_flags']
    shell:
        """
        mkdir -p plots/nanoplot
        (NanoPlot {params.flags} -o plots/nanoplot -p {wildcards.tag}_fastq_ --fastq_rich {input[0]} || NanoPlot {params.flags} -o plots/nanoplot -p {wildcards.tag}_fastq_ --fastq {input[0]})
        """

# can cause rule rerun for UMI plots
# rule NanoPlot_alignment_preConsensus:
#     input:
#         'sequences/UMI/{tag}_noConsensus.bam'
#     output:
#         'plots/nanoplot/{tag, [^\/_]*}_alignment_preConsensus_NanoStats.txt'
#     params:
#         flags = lambda wildcards: config['nanoplot_flags']
#     shell:
#         """
#         mkdir plots/nanoplot -p
#         NanoPlot {params.flags} -o plots/nanoplot -p {wildcards.tag}_alignment_preConsensus_ --bam {input[0]}
#         """

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

rule generate_barcode_ref:
    input:
        bam = lambda wildcards: expand('alignments/{tag}.bam', tag=config['barcode_fasta_to_tags'][wildcards.fasta]),
        barcode_info = lambda wildcards: config['runs'][config['barcode_fasta_to_tags'][wildcards.fasta][0]].get('barcode_info_csv', [])
    output:
        fasta = '{fasta}'
    params:
        references = lambda wildcards: [config['runs'][tag]['reference_aln'] for tag in config['barcode_fasta_to_tags'][wildcards.fasta]],
        tags = lambda wildcards: config['barcode_fasta_to_tags'][wildcards.fasta]
    wildcard_constraints:
        fasta = '|'.join([re.escape(f) for f in config.get('barcode_fasta_to_tags', {}).keys()]) if config.get('barcode_fasta_to_tags') else '__no_barcode_fastas__'
    script:
        'utils/generate_barcode_ref.py'

checkpoint demultiplex:
    input:
        aln = 'alignments/{tag}.bam',
        bai = 'alignments/{tag}.bam.bai',
        barcode_info = lambda wildcards: config['runs'][wildcards.tag].get('barcode_info_csv', []),
        partition_barcode_groups = lambda wildcards: config['runs'][wildcards.tag].get('partition_barcode_groups_csv', []),
        label_barcode_groups = lambda wildcards: config['runs'][wildcards.tag].get('label_barcode_groups_csv', []),
        barcode_fastas = lambda wildcards: [info['fasta'] for info in config['runs'][wildcards.tag].get('barcode_info', {}).values()]
    output:
        flag = touch('demux/.{tag, [^\/_]*}_demultiplex.done'),
        stats = 'demux/{tag, [^\/_]*}_demux-stats.csv'
    params:
        references_csv = lambda wildcards: config['runs'][wildcards.tag]['reference_csv'],
        threshold = lambda wildcards: config['runs'][wildcards.tag].get('demux_threshold', config.get('demux_threshold', 0)),
        screen_failures = lambda wildcards: config['runs'][wildcards.tag].get('demux_screen_failures', config.get('demux_screen_failures', True)),
        screen_no_group = lambda wildcards: config['runs'][wildcards.tag].get('demux_screen_no_group', config.get('demux_screen_no_group', True)),
        max_references = lambda wildcards: config.get('demux_max_references', False)
    script:
        'utils/demux.py'

rule merge_demux_stats:
    input:
        expand('demux/{tag}_demux-stats.csv', tag=[tag for tag in config['runs']])
    output:
        'demux-stats.csv'
    run:
        import pandas as pd
        dfs = [pd.read_csv(f) for f in input]
        combined = pd.concat(dfs)
        combined.to_csv(output[0], index=False)

# Proxy rule: declares demux BAM files as outputs after checkpoint completes
# The demultiplex checkpoint creates BAM files as side effects; this rule
# allows Snakemake to include them in the DAG, which is necessary for enrichment analysis.
rule demux_bam_proxy:
    input:
        flag = ancient('demux/.{tag}_demultiplex.done') # ancient prevents premature deletion
    output:
        bam = 'demux/{tag}_{barcodes}.bam'
    run:
        pass

rule index_demuxed:
    input:
        'demux/{tag}_{barcodes}.bam'
    output:
        'demux/{tag, [^\/]*}_{barcodes}.bam.bai'
    shell:
        """
        samtools index {input}
        """

rule plot_demux:
    input:
        CSV = 'demux/{tag}_demux-stats.csv'
    output:
        plot = 'plots/{tag, [^\/]*}_demux.html'
    params:
        barcodeInfo = lambda wildcards: config['runs'][wildcards.tag]['barcode_info'],
        barcodeGroups = lambda wildcards: config['runs'][wildcards.tag].get('barcodeGroups', {}),
        distribution_split_by_reference = lambda wildcards: config.get('distribution_split_by_reference', False)
    script:
        'utils/plot_demux.py'

def genotype_enrichment_input(wildcards):
    """
    Get genotype file inputs for genotype enrichment calculation.
    Finds all tag/barcode combinations used in this timepoint sample's series.
    Returns list of genotype CSV paths.

    This function is checkpoint-aware: it uses get_demuxed_barcodes_timepoint()
    to wait for demultiplex checkpoints and discover which barcodes actually exist.
    """
    timepoint_sample = wildcards.timepoint_sample
    if timepoint_sample not in config['timepointsInfo']:
        return []

    tag_barcode_tp = config['timepointsInfo'][timepoint_sample]['tag_barcode_tp']

    # Use get_demuxed_barcodes_timepoint to get only the barcodes that actually exist after demux
    # This function waits for checkpoints and filters out any barcodes that weren't produced
    available_tag_barcodes = get_demuxed_barcodes_timepoint(tag_barcode_tp.keys())

    genotype_files = []
    for (tag, barcode) in available_tag_barcodes:
        genotype_file = f'mutation_data/{tag}/{barcode}/{tag}_{barcode}_genotypes.csv'
        genotype_files.append(genotype_file)

    return genotype_files

# Prepare processed timepoints CSV for demux enrichment
# Adds -rep suffixes to duplicate sample names, keeps tag_barcodeGroup identifiers
rule prepare_timepoints_demux_csv:
    input:
        timepoints = lambda wildcards: config['timepoints']
    output:
        processed = 'enrichment/.timepoints_demux.csv'
    run:
        # Process timepoints CSV to add group/replicate columns and -rep suffixes
        processed_df = process_timepoints_csv(input.timepoints)

        # Save processed file
        processed_df.to_csv(output.processed, index=False)

# Prepare processed timepoints CSV for genotype enrichment
# Adds -rep suffixes and replaces tag_barcodeGroup with genotype file paths
rule prepare_timepoints_genotype_csv:
    input:
        timepoints = lambda wildcards: config['timepoints'],
        genotypes = lambda wildcards: [
            f'mutation_data/{tag}/{barcode}/{tag}_{barcode}_genotypes.csv'
            for sample in config.get('timepointsInfo', {}).keys()
            for (tag, barcode) in config['timepointsInfo'][sample].get('tag_barcode_tp', {}).keys()
        ] if 'timepointsInfo' in config else []
    output:
        processed = 'enrichment/.timepoints_genotype.csv'
    run:
        # Process timepoints CSV to add group/replicate columns and -rep suffixes
        processed_df = process_timepoints_csv(input.timepoints)

        # Replace tag_barcodeGroup identifiers with genotype file paths
        # Iterate through timepoint columns (skip group, replicate, sample_label)
        for col in processed_df.columns[3:]:  # Timepoint columns start at index 3
            processed_df[col] = processed_df[col].apply(
                lambda x: f'mutation_data/{x.split("_")[0]}/{x.split("_")[1]}/{x}_genotypes.csv'
                if pd.notna(x) and isinstance(x, str) and '_' in x else x
            )

        # Save processed file
        processed_df.to_csv(output.processed, index=False)

if config.get('enrichment_type') == 'demux':

    # Calculate demux enrichment scores for a timepoint sample
    # Demux enrichment uses barcode counts from demux-stats.csv
    rule demux_enrichment:
        input:
            CSV = 'demux-stats.csv',
            timepoints = 'enrichment/.timepoints_demux.csv'
        output:
            scores = 'enrichment/{timepoint_sample, [^\/]*}_demux-enrichment-scores.csv'
        params:
            mode = 'demux',
            timepoint_sample = lambda wildcards: wildcards.timepoint_sample,
            screen_failures = lambda wildcards: config.get('demux_screen_failures', False),
            reference_entity = lambda wildcards: config.get('enrichment_reference', '')
        threads: max(workflow.cores-1,1)
        script:
            'utils/enrichment.py'

elif config.get('enrichment_type') == 'genotype':

    # Calculate genotype enrichment scores for a timepoint sample
    # Genotype enrichment uses mutation data from genotype CSV files
    rule genotype_enrichment:
        input:
            timepoints = 'enrichment/.timepoints_genotype.csv',
            genotypes = genotype_enrichment_input  # Track genotype file dependencies
        output:
            scores = 'enrichment/{timepoint_sample, [^\/]*}_genotype-enrichment-scores.csv'
        params:
            mode = 'genotype',
            timepoint_sample = lambda wildcards: wildcards.timepoint_sample,
            reference_entity = lambda wildcards: config.get('enrichment_reference', '')
        threads: max(workflow.cores-1,1)
        script:
            'utils/enrichment.py'

def get_samples_for_group(group):
    """
    Get all sample_labels that belong to a given group.

    Args:
        group (str): The group name (e.g., 'BLF1531')

    Returns:
        list: Sample labels belonging to this group (e.g., ['BLF1531-rep1', 'BLF1531-rep2'])
              or ['BLF1531'] if no replicates
    """
    matching_samples = []
    if 'timepointsInfo' in config:
        for sample in config['timepointsInfo'].keys():
            sample_group = config['timepointsInfo'][sample].get('group', '')
            if sample_group == group:
                matching_samples.append(sample)
    return matching_samples


def enrichment_mean_input(wildcards):
    """
    Determine input enrichment scores files based on output mean file path.

    Handles replicate aggregation:
    - Input: Individual replicate score files (e.g., TrpB-IE5-rep1_*-enrichment-scores.csv, TrpB-IE5-rep2_*-enrichment-scores.csv)
    - Output: Aggregated mean file (e.g., TrpB-IE5_*-enrichment-scores-mean.csv, without -rep suffix)

    If only one sample exists (no replicates), uses that single score file directly.
    """
    mode = 'demux' if 'demux-enrichment' in wildcards.enrichment_path else 'genotype'

    # Extract base name from path (e.g., "TrpB-IE5_demux-enrichment-scores-mean" -> "TrpB-IE5")
    timepoint_base = wildcards.enrichment_path.split(f'_{mode}-enrichment')[0]

    # Find all timepoint samples that match this group
    matching_samples = get_samples_for_group(timepoint_base)

    if not matching_samples:
        raise ValueError(f"No timepoint samples found matching base name '{timepoint_base}' in timepointsInfo. "
                        f"Expected to find '{timepoint_base}' or '{timepoint_base}-rep1', etc.")

    # Return list of score files for all matching samples
    return [f'enrichment/{sample}_{mode}-enrichment-scores.csv' for sample in matching_samples]


def genotype_enrichment_scores_input(wildcards):
    group = wildcards.timepointsGroup
    matching_samples = get_samples_for_group(group)
    if not matching_samples:
        raise ValueError(f"No timepoint samples found for group '{group}' in timepointsInfo.")
    return [f'enrichment/{sample}_genotype-enrichment-scores.csv' for sample in matching_samples]

# Unified rule for filtering enrichment scores and computing means across replicates
# Input: One or more replicate enrichment score files (e.g., TrpB-IE5-rep1_*-scores.csv, TrpB-IE5-rep2_*-scores.csv)
# Output: Single mean file aggregating replicates (e.g., TrpB-IE5_*-scores-mean.csv)
rule enrichment_mean:
    input:
        scores = enrichment_mean_input
    output:
        mean = 'enrichment/{enrichment_path, .+enrichment-scores-mean}.csv'
    params:
        SE_filter = lambda wildcards: config.get('enrichment_SE_filter', 0),
        t0_filter = lambda wildcards: config.get('enrichment_t0_filter', 0),
        score_filter = lambda wildcards: config.get('enrichment_score_filter', False),
        include_mean_normalized = lambda wildcards: 'genotype-enrichment' in wildcards.enrichment_path
    run:
        from utils.enrichment import enrichment_mean_filter
        # Concatenate all replicate score files
        import pandas as pd
        all_scores = []
        for score_file in input.scores:
            df = pd.read_csv(score_file, index_col=False)
            all_scores.append(df)
        scores_df = pd.concat(all_scores, ignore_index=True)

        # Calculate mean and apply filters
        enrichment_mean_filter(scores_df, SE_filter=params.SE_filter, t0_filter=params.t0_filter,
                             score_filter=params.score_filter, mean_csv=output.mean,
                             include_mean_normalized=params.include_mean_normalized)

# Plot enrichment scores
rule plot_enrichment:
    input:
        mean = 'enrichment/{timepoint_sample}_{mode}-enrichment-scores-mean.csv'
    output:
        plot = 'plots/{timepoint_sample, [^\/]*}_{mode, (demux|genotype)}-enrichment-scores.html'
    run:
        from utils.enrichment import plot_enrichment
        mean_df = pd.read_csv(input.mean, index_col=False)
        plot_enrichment(mean_df, output.plot)

# NOTE: in previous Snakemake versions a workaround that allowed this rule to have a dynamic number of outputs depending on
#       if AA analysis is performed was used. This is no longer possible so when AA analysis is not performed, empty AA outputs
#       will be generated
rule mutation_analysis:
    input:
        bam = lambda wildcards: expand('demux/{tag}_{{barcodes}}.bam', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam',
        bai = lambda wildcards: expand('demux/{tag}_{{barcodes}}.bam.bai', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam.bai'
    output:
        expand('mutation_data/{{tag, [^\/_]*}}/{{barcodes, [^\/_]*}}/{{tag}}_{{barcodes}}_{datatype}', datatype = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-mutations-aggregated.csv', 'NT-mutation-distribution.csv', 'AA-mutations-aggregated.csv', 'AA-mutation-distribution.csv'])
    params:
        NT_muts_of_interest = lambda wildcards: config['runs'][wildcards.tag].get('NT_muts_of_interest',''),
        AA_muts_of_interest = lambda wildcards: config['runs'][wildcards.tag].get('AA_muts_of_interest',''),
        analyze_seqs_with_indels = lambda wildcards: config.get('analyze_seqs_with_indels', True),
        mutations_frequencies_raw = lambda wildcards: config.get('mutations_frequencies_raw', False),
        quality_score_minimum = lambda wildcards: config.get('mutation_analysis_quality_score_minimum', 5)
    script:
        'utils/mutation_analysis.py'

def get_demuxed_barcodes(tag, bcGroupsDict):
    """
    tag:            string, run tag
    bcGroupsDict:   dictionary, barcodeGroups dict from config file for a specific tag,
                        if it exists, otherwise an empty dict if no sorting is needed

    grabs all barcodes for a given tag that were properly demultiplexed
            if that tag was demultiplexed, then sorts according to the
            barcode groups in the config file"""

    if config['do_demux'][tag]:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split('demultiplex')[0]
        checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
        BCs = glob_wildcards(checkpoint_demux_files).BCs
        if not BCs:
            error_msg = f"""
{'='*80}
[ERROR] No demultiplexed barcode files found for tag `{tag}`
Expected files matching: {checkpoint_demux_files}

This usually means demultiplexing failed to produce any output files or no reads
matched the expected barcode contexts.
Check demux-stats.csv for details.
{'='*80}
"""
            raise RuntimeError(error_msg)
        out = sort_barcodes(BCs, bcGroupsDict)
    else:
        out = ['all']
    return out

def get_demuxed_barcodes_timepoint(tag_bcs):
    """
    filters a list of (tag,bc) tuples for availability after demux

    tag_bcs ( list(tuples) ): (tag, bc) tuples that will be checked for in the demux outputs.
        If any are not within the demux outputs, a warning will be pushed but the pipeline will proceed
    """
    out = []
    allowedBCsDict = {}
    for tag, bc in tag_bcs:
        if tag not in allowedBCsDict:
            allowedBCsDict[tag] = get_demuxed_barcodes(tag, {})
        if bc in allowedBCsDict[tag]: # add (tag,bc) tuple to output if that tag,bc combination is output by demux
            out.append( (tag,bc) )
        else:
            print(f'[NOTICE] tag_bc combination "{tag}_{bc}" was not demultiplexed. Not including in timepoint outputs')
    return out

def mut_stats_input(wildcards):
    datatypes = ['genotypes.csv', 'failures.csv', 'NT-mutations-aggregated.csv', 'NT-mutation-distribution.csv']
    if config['do_AA_mutation_analysis'][wildcards.tag]:
        datatypes.extend(['AA-mutations-aggregated.csv', 'AA-mutation-distribution.csv'])
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
	params:
		do_aa_analysis = lambda wildcards: config['do_AA_mutation_analysis'][wildcards.tag]
	script:
		'utils/mutation_statistics.py'

# allows for tolerance of tags that do not pass the demux checkpoint
def merge_mut_stats_input(wildcards):
    out = []

    for tag in config['runs']:
        if config['do_NT_mutation_analysis'][tag]:
            if config['do_demux'][tag]:
                checkpoint_demux_output = checkpoints.demultiplex.get(tag=tag).output[0]
                checkpoint_demux_prefix = checkpoint_demux_output.split('demultiplex')[0]
                checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
                if len(glob_wildcards(checkpoint_demux_files).BCs) > 0:
                    out.append(f'mutation_data/{tag}/{tag}_mutation-stats.csv')
            else:
                out.append(f'mutation_data/{tag}/{tag}_mutation-stats.csv')
    return out

rule merge_mut_stats:
    input:
        merge_mut_stats_input
    output:
        'mutation-stats.csv'
    run:
        import pandas as pd
        dfs = [pd.read_csv(f) for f in input]
        combined = pd.concat(dfs)
        combined.to_csv(output[0], index=False)

def get_timepoint_plots_all_input(wildcards):
    out = []
    # if not config.get('timepoints', {}):
    #     print('[WARNING] All timepoints rule invoked but timepoints not specified in config file. No timepoint plots will be generated.')
    # else:
    # other timepoints outputs are specific to each row in the timepoints file
    timepoint_rows = [row for row in config['timepointsInfo'] if len(get_demuxed_barcodes_timepoint( config['timepointsInfo'][row]['tag_barcode_tp'].keys() )) > 0] # only include timepoints that have at least one demuxed barcode
    print_flag = False
    plot_types = ['mutation-distribution-violin', 'genotypes-distribution']
    if config.get('plot_mutations-aggregated', False):
        plot_types.append('AA-mutations-aggregated')
    if config.get('genotypes2D_plot_groups', False):
        plot_types.append('genotypes2D')
    out.extend(expand('plots/timepoints/{timepointSample}_{plot_type}.html', timepointSample=timepoint_rows, plot_type=plot_types))    # each of these outputs includes data for a single row in the timepoints file
    timepointSample_NTorAA = []
    tags = []
    for timepointSample in timepoint_rows:
        tag = config['timepointsInfo'][timepointSample]['tag']
        if tag not in tags:
            tags.append(tag)
        timepointSample_NTorAA.extend(expand('{TS}_{NTorAA}', TS=timepointSample, NTorAA= ['NT','AA'] if config['do_AA_mutation_analysis'][ tag ] else ['NT']))
    # out.extend(expand('plots/{tag}_mutation-rates-mut-grouped.html', tag=[t for t in tags if config['do_NT_mutation_analysis'][t]]))                 # mutation rates includes data for all rows in the timepoints file
    dist_types = ['mutation']
    if config['plot_hamming-distance-distribution']:
        dist_types.append('hamming-distance')
    out.extend(expand('plots/timepoints/{timepointSample_NTorAA}-{distType}-distribution.html', timepointSample_NTorAA=timepointSample_NTorAA, distType=dist_types))

    return out

rule timepoint_plots_all:
    input:
        get_timepoint_plots_all_input
    output:
        'plots/.all_timepoints.done'
    shell:
        'touch {output}'

rule genotypes_distribution:
    input:
        genotypes = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv'
    output:
        distribution = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes-distribution.csv'
    run:
        import pandas as pd
        import numpy as np

        df = pd.read_csv(input.genotypes)
        max_value = df['count'].max()
        bins = np.concatenate(([0], 10 ** np.arange(0, np.ceil(np.log10(max_value)) + 1))) # manually add in 0 bin for unique genotypes
        binned_sums = df['count'].groupby(pd.cut(df['count'], bins=bins, right=True)).sum()
        df = pd.DataFrame({'bin': binned_sums.index, 'total reads': binned_sums.values})
        df['log(max of bin)'] = df['bin'].apply(lambda x: int(np.log10(x.right)))
        df['proportion of all reads'] = df['total reads'] / df['total reads'].sum()
        df['cumulative proportion of all reads'] = df['proportion of all reads'].cumsum()
        df = df[['log(max of bin)', 'total reads', 'proportion of all reads', 'cumulative proportion of all reads']]
        df.to_csv(output.distribution, index=False)


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
        timepoints = lambda wildcards: config['timepoints']
    output:
        boxplot_mut_grouped = 'plots/{tag, [^\/]*}_mutation-rates-mut-grouped.html',
        boxplot_plot_sample_grouped = 'plots/{tag, [^\/]*}_mutation-rates-sample-grouped.html',
        heatmap = 'plots/{tag, [^\/]*}_mutation-rates-heatmap.html',
        CSV_all_rates = 'mutation_data/{tag, [^\/]*}/{tag}_mutation-rates.csv',
        CSV_summary = 'mutation_data/{tag, [^\/]*}/{tag}_mutation-rates-summary.csv',
    params:
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r')
    script:
        'utils/plot_mutation_rate.py'

rule plot_mutations_aggregated_tag:
    input:
        genotypes = lambda wildcards: expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{NTorAA}-mutations-aggregated.csv',
                                                    tag=wildcards.tag,
                                                    barcodes=get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {})),
                                                    NTorAA = wildcards.NTorAA ),
        mut_stats = 'mutation_data/{tag}/{tag}_mutation-stats.csv'
    output:
        all_muts = 'plots/{tag, [^\/_]*}_{NTorAA, [^\/_]*}-mutations-aggregated.html',
        most_frequent = 'plots/{tag, [^\/_]*}_{NTorAA, [^\/_]*}-mutations-aggregated-common.html',
        least_frequent = 'plots/{tag, [^\/_]*}_{NTorAA, [^\/_]*}-mutations-aggregated-rare.html'
    params:
        mutations_frequencies_raw = lambda wildcards: config.get('mutations_frequencies_raw', False),
        number_of_positions = lambda wildcards: config.get('mutations_frequencies_number_of_positions', 20),
        heatmap = lambda wildcards: config.get('mutations_frequencies_heatmap', False),
        labels = lambda wildcards, input: [in_file.split('/')[-2] for in_file in input],
        title = lambda wildcards: wildcards.tag,
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        split_by_reference = lambda wildcards: config.get('mutations_aggregated_split_by_reference', 10)
    script:
        'utils/plot_mutations_frequencies.py'

rule plot_mutations_aggregated_timepointSample:
    input:
        genotypes = lambda wildcards: expand('mutation_data/{tag_bc_tag_bc}_{NTorAA}-mutations-aggregated.csv',
                                                    tag_bc_tag_bc = [f'{tag}/{barcodes}/{tag}_{barcodes}' for tag,barcodes in get_demuxed_barcodes_timepoint( config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].keys() )],
                                                    NTorAA = wildcards.NTorAA ),
        mut_stats = 'mutation-stats.csv'
    output:
        all_muts = 'plots/timepoints/{timepointsSample, [^\/_]*}_{NTorAA, [^\/_-]*}-mutations-aggregated.html',
        most_frequent = 'plots/timepoints/{timepointsSample, [^\/_]*}_{NTorAA, [^\/_-]*}-mutations-aggregated-common.html',
        least_frequent = 'plots/timepoints/{timepointsSample, [^\/_]*}_{NTorAA, [^\/_-]*}-mutations-aggregated-rare.html'
    params:
        mutations_frequencies_raw = lambda wildcards: config.get('mutations_frequencies_raw', False),
        number_of_positions = lambda wildcards: config.get('mutations_frequencies_number_of_positions', 20),
        heatmap = lambda wildcards: config.get('mutations_frequencies_heatmap', False),
        labels = lambda wildcards: list(config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].values()),
        title = lambda wildcards: wildcards.timepointsSample,
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        split_by_reference = lambda wildcards: config.get('mutations_aggregated_split_by_reference', 10)
    script:
        'utils/plot_mutations_frequencies.py'

rule plot_mutations_aggregated:
    input:
        frequencies = 'mutation_data/{tag}_{barcodes}_{NTorAA}-mutations-aggregated.csv',
        mutStats = 'mutation_data/{tag}/{tag}_mutation-stats.csv'
    output:
        all_muts = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA, [^\/_]*}-mutations-aggregated.html',
        most_frequent = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA, [^\/_]*}-mutations-aggregated-common.html',
        least_frequent = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA, [^\/_]*}-mutations-aggregated-rare.html',
    params:
        mutations_frequencies_raw = lambda wildcards: config.get('mutations_frequencies_raw', False),
        number_of_positions = lambda wildcards: config.get('mutations_frequencies_number_of_positions', 20),
        heatmap = lambda wildcards: config.get('mutations_frequencies_heatmap', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r'),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        split_by_reference = lambda wildcards: config.get('mutations_aggregated_split_by_reference', 10)
    script:
        'utils/plot_mutations_frequencies.py'

rule hamming_distance:
    input:
        genotypesCSV = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv',
        ref_csv = lambda wildcards: config['runs'][wildcards.tag]['reference_csv']
    output:
        # HDmatrixCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA}-hamming-distance-matrix.csv', # not currently using and takes up a lot of disk space
        HDdistCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA}-hamming-distance-distribution.csv'
    params:
        downsample = lambda wildcards: config.get('hamming_distance_distribution_downsample', False),
        analyze_seqs_with_indels = lambda wildcards: config.get('analyze_seqs_with_indels', True)
    run:
        import pandas as pd
        from utils.common import dist_to_DF, load_references_from_csv
        from utils.SequenceAnalyzer import SequenceAnalyzer

        # Load references
        references_dict, errors, _, _, _ = load_references_from_csv(input.ref_csv)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            raise ValueError("Failed to load references from CSV")

        # Load genotypes to determine which references are present
        genotypes_df = pd.read_csv(input.genotypesCSV)

        # Process each reference separately
        hd_dist_dfs = []
        for ref_name in genotypes_df['reference_name'].unique():
            # Filter genotypes for this reference
            ref_genotypes = genotypes_df[genotypes_df['reference_name'] == ref_name]

            # Get reference sequences for this reference [alignment_seq, NT_seq, AA_seq]
            ref_data = references_dict[ref_name]
            reference_sequences = [ref_data['alignment_seq'], ref_data['NT_seq'], ref_data['coding_seq']]

            # Create SequenceAnalyzer with filtered genotypes
            data = SequenceAnalyzer(reference_sequences=reference_sequences,
                                   genotypes_df=ref_genotypes,
                                   exclude_indels=(not params.analyze_seqs_with_indels))

            # Calculate hamming distance
            matrixDF, HDdistDF = data.HD_matrix_and_dist(NTorAA=wildcards.NTorAA, downsample=params.downsample)

            # Convert to DataFrame and add reference_name column
            HD_dist_DF = dist_to_DF(HDdistDF, f'{wildcards.NTorAA} hamming distance', 'sequence pairs')
            HD_dist_DF['reference_name'] = ref_name

            hd_dist_dfs.append(HD_dist_DF)

        # Concatenate all references and save
        combined_hd_dist = pd.concat(hd_dist_dfs, ignore_index=True)
        combined_hd_dist.to_csv(output.HDdistCSV, index=False)

def plot_mutations_distribution_input(wildcards):
    if config['do_demux'][wildcards.tag]:
        checkpoint_demux_output = checkpoints.demultiplex.get(tag=wildcards.tag).output[0]
        checkpoint_demux_prefix = checkpoint_demux_output.split(f'demultiplex')[0]
        checkpoint_demux_files = checkpoint_demux_prefix.replace('.','') + '{BCs}.bam'
        return expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{NTorAA}-hamming-distance-distribution.csv', tag=wildcards.tag, barcodes=glob_wildcards(checkpoint_demux_files).BCs, NTorAA = wildcards.NTorAA)
    else:
        return expand('mutation_data/{tag}/all/{tag}_all_{NTorAA}-hamming-distance-distribution.csv', tag=wildcards.tag, NTorAA=wildcards.NTorAA)

rule plot_distribution:
    input:
        'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{NTorAA}-{distType}-distribution.csv'
    output:
        plot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA, [^\/_]*}-{distType, [^\/_]*}-distribution.html'
    params:
        labels = lambda wildcards: [f"{wildcards.tag}, {wildcards.barcodes}"],
        title = lambda wildcards: wildcards.tag,
        legend_label = False,
        background = False,
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        x_max = lambda wildcards: config.get('distribution_x_range', False),
        y_max = lambda wildcards: config.get('distribution_y_range', False),
        split_by_reference = lambda wildcards: config.get('split_by_reference', False)
    script:
        'utils/plot_distribution.py'

rule plot_distribution_tag:
    input:
        lambda wildcards: expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{NTorAA}-{distType}-distribution.csv',
                                tag = wildcards.tag,
                                barcodes = get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {})),
                                NTorAA = wildcards.NTorAA,
                                distType = wildcards.distType)
    output:
        plot = 'plots/{tag, [^\/_]*}_{NTorAA, [^\/_]*}-{distType, [^\/_]*}-distribution.html'
    params:
        labels = lambda wildcards, input: [in_file.split('/')[-2] for in_file in input],
        title = lambda wildcards: wildcards.tag,
        legend_label = 'barcode group',
        background = lambda wildcards: config.get('background', False),
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        x_max = lambda wildcards: config.get('distribution_x_range', False),
        y_max = lambda wildcards: config.get('distribution_y_range', False),
        split_by_reference = lambda wildcards: config.get('split_by_reference', False)
    script:
        'utils/plot_distribution.py'

rule plot_distribution_timepointSample:
    input:
        lambda wildcards: expand('mutation_data/{tag_bc_tag_bc}_{NTorAA}-{distType}-distribution.csv',
                                tag_bc_tag_bc = [f'{tag}/{barcodes}/{tag}_{barcodes}' for tag,barcodes in get_demuxed_barcodes_timepoint( config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].keys() )],
                                NTorAA = wildcards.NTorAA,
                                distType = wildcards.distType)
    output:
        plot = 'plots/timepoints/{timepointsSample, [^\/_]*}_{NTorAA, [^\/_-]*}-{distType, [^\/_]*}-distribution.html'
    params:
        labels = lambda wildcards: list(config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].values()),
        title = lambda wildcards: wildcards.timepointsSample,
        legend_label = lambda wildcards: config['timepointsInfo'][wildcards.timepointsSample]['units'],
        background = lambda wildcards: config.get('background', False),
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        x_max = lambda wildcards: config.get('distribution_x_range', False),
        y_max = lambda wildcards: config.get('distribution_y_range', False),
        split_by_reference = lambda wildcards: config.get('split_by_reference', False)
    script:
        'utils/plot_distribution.py'

rule plot_genotypes_distribution_timepointSample:
    input:
        lambda wildcards: expand('mutation_data/{tag_bc_tag_bc}_genotypes-distribution.csv',
                                tag_bc_tag_bc = [f'{tag}/{barcodes}/{tag}_{barcodes}' for tag,barcodes in get_demuxed_barcodes_timepoint( config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].keys() )])
    output:
        plot = 'plots/timepoints/{timepointsSample, [^\/_]*}_genotypes-distribution.html'
    params:
        labels = lambda wildcards: list(config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].values()),
        title = lambda wildcards: wildcards.timepointsSample,
        legend_label = lambda wildcards: config['timepointsInfo'][wildcards.timepointsSample]['units'],
        background = lambda wildcards: config.get('background', False),
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        colormap = lambda wildcards: config.get('colormap', 'kbc_r'),
        x_max = lambda wildcards: config.get('distribution_x_range', False),
        y_max = lambda wildcards: config.get('distribution_y_range', False),
        split_by_reference = lambda wildcards: config.get('split_by_reference', False)
    script:
        'utils/plot_distribution.py'

rule plot_violin_distribution_tag:
    input:
        'mutation_data/{tag}/{tag}_genotypes.csv'
    output:
        'plots/{tag, [^\/_]*}_mutation-distribution-violin.html'
    params:
        group_col = 'barcode_group',
        x_label = 'sample',
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r'),
        background = lambda wildcards: config.get('background', False)
    script:
        'utils/plot_mutation_violin_distribution.py'

rule plot_violin_distribution_timepointSample:
    input:
        'mutation_data/timepoints/{timepointsSample}_merged-timepoint_genotypes.csv'
    output:
        'plots/timepoints/{timepointsSample, [^\/_]*}_mutation-distribution-violin.html'
    params:
        group_col = 'timepoint',
        x_label = lambda wildcards: config['timepointsInfo'][wildcards.timepointsSample]['units'],
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r'),
        background = False
    script:
        'utils/plot_mutation_violin_distribution.py'

rule reduce_genotypes_dimensions:
    input:
        genotypes = '{dir}/{tag}_{barcodes}_genotypes.csv'
    output:
        reduced = '{dir}/{tag, [^\/_]*}_{barcodes, [^\/]*}_genotypes-reduced-dimensions.csv'
    params:
        ref_csv = lambda wildcards: config['runs'][wildcards.tag].get('reference_csv', '')
    run:
        import pandas as pd
        from utils.SequenceAnalyzer import SequenceAnalyzer
        from utils.common import load_references_from_csv

        # Load genotypes
        genotypes_df = pd.read_csv(input.genotypes)

        # Add reference_name column if not present (backward compatibility)
        if 'reference_name' not in genotypes_df.columns:
            genotypes_df['reference_name'] = 'N/A'

        # Load references
        references_dict, errors, _, _, _ = load_references_from_csv(params.ref_csv)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            raise ValueError("Failed to load references from CSV")

        # Process each reference separately
        reduced_dfs = []
        for ref_name in genotypes_df['reference_name'].unique():
            # Filter genotypes for this reference
            ref_genotypes = genotypes_df[genotypes_df['reference_name'] == ref_name]

            # Get reference sequences [alignment_seq, NT_seq, coding_seq]
            ref_data = references_dict[ref_name]
            reference_sequences = [ref_data['alignment_seq'], ref_data['NT_seq'], ref_data['coding_seq']]

            # Create SequenceAnalyzer with filtered genotypes
            sequences = SequenceAnalyzer(reference_sequences=reference_sequences, genotypes_df=ref_genotypes, exclude_indels=False)
            sequences.assign_dimension_reduction('NT')
            if sequences.do_AA_analysis:
                sequences.assign_dimension_reduction('AA')

            reduced_dfs.append(sequences.genotypes)

        # Concatenate all references
        combined = pd.concat(reduced_dfs, ignore_index=True)
        combined.to_csv(output.reduced, index=False)

rule reduce_timepoint_genotypes_dimensions:
    input:
        genotypes = 'mutation_data/timepoints/{timepointsSample}_merged-timepoint_genotypes.csv'
    output:
        reduced = 'mutation_data/timepoints/{timepointsSample, [^\/_]*}_merged-timepoint_genotypes-reduced-dimensions.csv'
    params:
        ref_csv = lambda wildcards: config['timepointsInfo'][wildcards.timepointsSample].get('reference_csv', '')
    run:
        import pandas as pd
        from utils.SequenceAnalyzer import SequenceAnalyzer
        from utils.common import load_references_from_csv

        # Load genotypes
        genotypes_df = pd.read_csv(input.genotypes)

        # Add reference_name column if not present (backward compatibility)
        if 'reference_name' not in genotypes_df.columns:
            genotypes_df['reference_name'] = 'N/A'

        # Load references
        references_dict, errors, _, _, _ = load_references_from_csv(params.ref_csv)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            raise ValueError("Failed to load references from CSV")

        # Process each reference separately
        reduced_dfs = []
        for ref_name in genotypes_df['reference_name'].unique():
            # Filter genotypes for this reference
            ref_genotypes = genotypes_df[genotypes_df['reference_name'] == ref_name]

            # Get reference sequences [alignment_seq, NT_seq, coding_seq]
            ref_data = references_dict[ref_name]
            reference_sequences = [ref_data['alignment_seq'], ref_data['NT_seq'], ref_data['coding_seq']]

            # Create SequenceAnalyzer with filtered genotypes
            sequences = SequenceAnalyzer(reference_sequences=reference_sequences, genotypes_df=ref_genotypes, exclude_indels=False)
            sequences.assign_dimension_reduction('NT')
            if sequences.do_AA_analysis:
                sequences.assign_dimension_reduction('AA')

            reduced_dfs.append(sequences.genotypes)

        # Concatenate all references
        combined = pd.concat(reduced_dfs, ignore_index=True)
        combined.to_csv(output.reduced, index=False)

rule merge_tag_genotypes:
    input:
        genotypeCSVs = lambda wildcards: expand( 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv', tag=wildcards.tag, barcodes=get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {})) ),
    output:
        mergedGenotypes = temp('mutation_data/{tag}/{tag}_genotypes.csv')
    run:
        import pandas as pd
        DFs = []
        for csv in input.genotypeCSVs:
            df = pd.read_csv(csv)
            barcode_group = csv.split('_')[-2]
            df['barcode_group'] = barcode_group
            DFs.append(df)
        pd.concat(DFs).to_csv(output.mergedGenotypes, index=False)

def get_demux_seq_ids_path(wildcards):
    """Get seq-IDs path for demux enrichment merge - uses first sample's first barcode."""
    sample = get_samples_for_group(wildcards.timepointsGroup)[0]
    tag_barcode_tp = config['timepointsInfo'][sample]['tag_barcode_tp']
    tag, barcode = list(tag_barcode_tp.keys())[0]
    return f'mutation_data/{tag}/{barcode}/{tag}_{barcode}_seq-IDs.csv'

rule merge_timepoint_demux_enrichment_with_genotypes:
    input:
        sample_genotypes = lambda wildcards: [f'mutation_data/timepoints/{s}_merged-timepoint_genotypes-reduced-dimensions.csv'
                                               for s in get_samples_for_group(wildcards.timepointsGroup)],
        seq_ids = get_demux_seq_ids_path,
        enrichment = 'enrichment/{timepointsGroup}_demux-enrichment-scores-mean.csv'
    output:
        genotypes_enrichment = 'mutation_data/timepoints/{timepointsGroup, [^\/_]*}_merged-group_genotypes-reduced-dimensions-demux-enrichment-mean.csv'
    params:
        filter_missing_replicates = lambda wildcards: config.get('enrichment_missing_replicates_filter', True),
        dedupe_columns = ['genotype_ID', 'reference_name']
    script:
        'utils/merge_enrichment.py'

rule reduce_genotype_enrichment_dimensions:
    input:
        mean = 'enrichment/{timepointsGroup}_genotype-enrichment-scores-mean.csv',
        scores = genotype_enrichment_scores_input
    output:
        reduced = 'mutation_data/timepoints/{timepointsGroup, [^\/_]*}_merged-group_genotypes-reduced-dimensions-genotype-enrichment-mean.csv'
    params:
        ref_csv = lambda wildcards: config['timepointsInfo'][get_samples_for_group(wildcards.timepointsGroup)[0]].get('reference_csv', '')
    run:
        import pandas as pd
        import sys
        from utils.SequenceAnalyzer import SequenceAnalyzer
        from utils.common import load_references_from_csv

        mean_df = pd.read_csv(input.mean)

        # Scores files have mutation detail columns that mean file lacks
        all_scores = [pd.read_csv(f) for f in input.scores]
        scores_df = pd.concat(all_scores, ignore_index=True)

        mutation_cols = ['entity_id', 'reference_name', 'NT_substitutions', 'NT_insertions', 'NT_deletions',
                        'AA_substitutions_nonsynonymous', 'AA_substitutions_synonymous']
        mutation_cols = [col for col in mutation_cols if col in scores_df.columns]
        mutations_unique = scores_df[mutation_cols].drop_duplicates(subset=['entity_id'])

        genotypes_df = mean_df.merge(mutations_unique, on='entity_id', how='left')
        genotypes_df['count'] = 1  # Placeholder for SequenceAnalyzer

        references_dict, errors, _, _, _ = load_references_from_csv(params.ref_csv)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            raise ValueError("Failed to load references from CSV")

        reduced_dfs = []
        for ref_name in genotypes_df['reference_name'].unique():
            ref_genotypes = genotypes_df[genotypes_df['reference_name'] == ref_name]
            ref_data = references_dict[ref_name]
            reference_sequences = [ref_data['alignment_seq'], ref_data['NT_seq'], ref_data['coding_seq']]

            sequences = SequenceAnalyzer(reference_sequences=reference_sequences, genotypes_df=ref_genotypes, exclude_indels=False)
            sequences.assign_dimension_reduction('NT')
            if sequences.do_AA_analysis:
                sequences.assign_dimension_reduction('AA')

            reduced_dfs.append(sequences.genotypes)

        combined = pd.concat(reduced_dfs, ignore_index=True)
        combined.to_csv(output.reduced, index=False)

rule plot_genotypes2D:
    input:
        genotypesReduced = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes-reduced-dimensions.csv'
    output:
        genotypes2Dscatter = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_genotypes2D.html',
        genotypes2Dhexbins = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_genotypes2Dhexbins.html'
    params:
        downsample = lambda wildcards: config.get('genotypes2D_plot_downsample', False),
        plot_AA = lambda wildcards: config.get('genotypes2D_plot_AA', False) if config['do_AA_mutation_analysis'][wildcards.tag] else False,
        size_column = lambda wildcards: config.get('genotypes2D_plot_point_size_col', 'count'),
        size_range = lambda wildcards: config.get('genotypes2D_plot_point_size_range', '10, 30'),
        color_column = lambda wildcards: config.get('genotypes2D_plot_point_color_col', 'NT_substitutions_count'),
        background = lambda wildcards: config.get('background', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r')
    script:
        'utils/plot_genotypes_2d.py'

rule plot_genotypes2D_bcGroup:
    input:
        genotypesReduced = 'mutation_data/{tag}/{tag}_genotypes-reduced-dimensions.csv'
    output:
        genotypes2Dscatter = 'plots/{tag, [^\/_]*}_genotypes2D.html',
        genotypes2Dhexbins = 'plots/{tag, [^\/_]*}_genotypes2Dhexbins.html'
    params:
        downsample = lambda wildcards: config.get('genotypes2D_plot_downsample', False),
        plot_AA = lambda wildcards: config.get('genotypes2D_plot_AA', False) if config['do_AA_mutation_analysis'][wildcards.tag] else False,
        size_column = lambda wildcards: config.get('genotypes2D_plot_point_size_col', 'count'),
        size_range = lambda wildcards: config.get('genotypes2D_plot_point_size_range', '10, 30'),
        color_column = lambda wildcards: config.get('genotypes2D_plot_point_color_col', 'barcode_group'),
        background = lambda wildcards: config.get('background', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r')
    script:
        'utils/plot_genotypes_2d.py'

rule merge_timepoint_genotypes:
    input:
        genotypeCSVs = lambda wildcards: expand('mutation_data/{tag_barcodes}_genotypes.csv', tag_barcodes=[f'{tag}/{barcodes}/{tag}_{barcodes}' for tag,barcodes in get_demuxed_barcodes_timepoint( config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp'].keys() )]),
    output:
        mergedGenotypes = temp('mutation_data/timepoints/{timepointsSample, [^\/_]*}_merged-timepoint_genotypes.csv')
    params:
        tpInfo = lambda wildcards: config['timepointsInfo'][wildcards.timepointsSample]['tag_barcode_tp']
    run:
        import pandas as pd
        timepoints = list(params.tpInfo.values())   # timepoint values are an amount of some unit of time like X generations
        DFs = []
        if len(input.genotypeCSVs) > 0:
            for i, csv in enumerate(input.genotypeCSVs):
                df = pd.read_csv(csv)
                df['timepoint'] = timepoints[i]
                DFs.append(df)
            out_df = pd.concat(DFs)
        else:
            raise InputException('No genotypes found for timepoints sample: ' + wildcards.timepointsSample)
        out_df = out_df[['timepoint'] + [c for c in out_df.columns if c != 'timepoint']]
        out_df.to_csv(output.mergedGenotypes, index=False)

rule plot_genotypes2D_timepointSample:
    input:
        genotypesReduced = 'mutation_data/timepoints/{timepointsSample}_merged-timepoint_genotypes-reduced-dimensions.csv'
    output:
        genotypes2Dscatter = 'plots/timepoints/{timepointsSample, [^\/_]*}_genotypes2D.html',
        genotypes2Dhexbins = 'plots/timepoints/{timepointsSample, [^\/_]*}_genotypes2Dhexbins.html'
    params:
        downsample = lambda wildcards: config.get('genotypes2D_plot_downsample', False),
        plot_AA = lambda wildcards: config.get('genotypes2D_plot_AA', False) if config['do_AA_mutation_analysis'][ config['timepointsInfo'][wildcards.timepointsSample]['tag'] ] else False,
        size_column = lambda wildcards: config.get('genotypes2D_plot_point_size_col', 'count'),
        size_range = lambda wildcards: config.get('genotypes2D_plot_point_size_range', '10, 30'),
        color_column = lambda wildcards: config.get('genotypes2D_plot_point_color_col', 'timepoint'),
        background = lambda wildcards: config.get('background', False),
        export_SVG = lambda wildcards: config.get('export_SVG', False),
        cmap = lambda wildcards: config.get('colormap', 'kbc_r')
    script:
        'utils/plot_genotypes_2d.py'

def plot_genotypes_2D_all_input(wildcards):
    out = []
    for tag in config['runs']:
        if config['do_NT_mutation_analysis'][tag]:
            if config['do_demux'][tag]:
                BCs = get_demuxed_barcodes(tag, config['runs'][tag].get('barcodeGroups', {}))
            else:
                BCs = ['all']
            out.extend( expand('plots/{tag}/{barcodes}/{tag}_{barcodes}_genotypes2D.html', tag=tag, barcodes=BCs) ) 
    return out

# plot genotypes 2D for each sample individually
rule plot_genotypes_2D_all:
    input:
        plot_genotypes_2D_all_input
    output:
        touch('plots/.all_genotypes_2D.done')

rule structure_heatmap:
    input:
        agg_csv = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_AA-mutations-aggregated.csv',
        structure_file = lambda wildcards: config['runs'][wildcards.tag].get('structure_file', '')
    output:
        html = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_structure-heatmap.html'
    params:
        colormap = lambda wildcards: config.get('colormap', 'bmy_r'),
        width = lambda wildcards: config.get('structure_heatmap_width', 800),
        height = lambda wildcards: config.get('structure_heatmap_height', 600),
        label_step = lambda wildcards: config.get('structure_heatmap_label_step', 50),
        min_identity = lambda wildcards: config.get('structure_heatmap_min_identity', 0.95)
    script:
        'utils/structure_heatmap.py'

rule run_dashboard_single_input:
    input:
        unpack(lambda wildcards: dashboard_input(wildcards=wildcards, config=config, single=True))
    params:
        port = config.get('dashboard_port', 3365),
        basedir = workflow.basedir,
        exclude_indels = lambda wildcards: '' if config.get('analyze_seqs_with_indels', True) else '--exclude_indels'
    shell:
        """
        panel serve --show --port {params.port} {params.basedir}/rules/utils/genotypes_dashboard.py --args --genotypes={input.genotypes} --reference={input.refFasta} --structure_file={input.structure_file} {params.exclude_indels}
        """

rule create_dashboard_datasets_csv:
    output:
        csv = f"{config.get('metadata_folder', 'metadata')}/.dashboard_datasets.csv"
    params:
        do_enrichment = config.get('do_enrichment', False),
        enrichment_type = config.get('enrichment_type', 'genotype'),
        timepointsInfo = config.get('timepointsInfo', {}),
        analyze_seqs_with_indels = config.get('analyze_seqs_with_indels', True)
    run:
        import os
        df = get_dashboard_datasets_df()
        os.makedirs(os.path.dirname(output.csv), exist_ok=True)
        df.to_csv(output.csv, index=False)

rule run_dashboard:
    input:
        datasets_csv = lambda wildcards: dashboard_input(wildcards=wildcards, config=config, single=False)['datasets_csv'],
        genotypes = lambda wildcards: get_dashboard_datasets_df()['genotypes_csv'].tolist(),
        structure_files = lambda wildcards: dashboard_input(wildcards=wildcards, config=config, single=False)['structure_files']
    params:
        port = config.get('dashboard_port', 3365),
        basedir = workflow.basedir
    shell:
        """
        panel serve --show --port {params.port} {params.basedir}/rules/utils/genotypes_dashboard.py --args --datasets={input.datasets_csv}
        """


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
