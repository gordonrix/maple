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
from natsort import natsorted

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

rule plot_demux:
    input:
        CSV = 'demux/{tag}_demux-stats.csv'
    output:
        plot = 'plots/{tag, [^\/]*}_demux.html'
    params:
        barcodeInfo = lambda wildcards: config['runs'][wildcards.tag]['barcodeInfo'],
        barcodeGroups = lambda wildcards: config['runs'][wildcards.tag].get('barcodeGroups', {})
    script:
        'utils/plot_demux.py'

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
        expand('mutation_data/{{tag, [^\/_]*}}/{{barcodes, [^\/_]*}}/{{tag}}_{{barcodes}}_{datatype}', datatype = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-mutation-frequencies.csv', 'NT-mutation-distribution.csv'])
    params:
        NT_muts_of_interest = lambda wildcards: config['runs'][wildcards.tag].get('NT_muts_of_interest','')
    script:
        'utils/mutation_analysis.py'

rule mutation_analysis:
    input:
        bam = lambda wildcards: expand('demux/{tag}_{{barcodes}}.bam', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam',
        bai = lambda wildcards: expand('demux/{tag}_{{barcodes}}.bam.bai', tag=wildcards.tag) if config['do_demux'][wildcards.tag] else f'alignments/{wildcards.tag}.bam.bai'
    output:
        expand('mutation_data/{{tag, [^\/_]*}}/{{barcodes, [^\/_]*}}/{{tag}}_{{barcodes}}_{datatype}', datatype = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-mutation-frequencies.csv', 'NT-mutation-distribution.csv', 'AA-mutation-frequencies.csv', 'AA-mutation-distribution.csv'])
    params:
        NT_muts_of_interest = lambda wildcards: config['runs'][wildcards.tag].get('NT_muts_of_interest',''),
        AA_muts_of_interest = lambda wildcards: config['runs'][wildcards.tag].get('AA_muts_of_interest','')
    script:
        'utils/mutation_analysis.py'

def mut_stats_input(wildcards):
    datatypes = ['alignments.txt', 'genotypes.csv', 'seq-IDs.csv', 'failures.csv', 'NT-mutation-frequencies.csv', 'NT-mutation-distribution.csv']
    if config['do_AA_mutation_analysis'][wildcards.tag]: datatypes.extend(['AA-mutation-frequencies.csv', 'AA-mutation-distribution.csv'])
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

def sort_barcodes(bcList, bcGroupsDict):
    """
    bcList:         list of strings, barcode groups
    bcGroupsList:   dict, barcode groups dictionary for a specific tag from the config file,
                        or an empty dictionary

    returns the same list of barcode sorted first according to user provided barcode groups,
        then by natural sorting (see natsort package)
    """
    defined = [bc for bc in bcGroupsDict.keys() if bc in bcList]
    undefined = natsorted([bc for bc in bcList if bc not in bcGroupsDict.keys()])
    return defined+undefined

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
    

def dms_view_input(wildcards):
    out = []
    for tag in config['runs']:
        out.extend( expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_AA-mutation-frequencies.csv', tag=tag, barcodes=get_demuxed_barcodes(tag, config['runs'][tag].get('barcodeGroups', {}) )) )
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

rule plot_mutations_frequencies:
    input:
        frequencies = lambda wildcards: expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{AAorNT}-mutation-frequencies.csv',
                                                    tag=wildcards.tag,
                                                    barcodes=get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {})),
                                                    AAorNT = wildcards.AAorNT ),
        mutStats = 'mutation_data/{tag}/{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/_]*}_{AAorNT, [^\/_]*}-mutation-frequencies.html'
    script:
        'utils/plot_mutations_frequencies.py'

rule plot_mutations_frequencies_barcodeGroup:
    input:
        frequencies = 'mutation_data/{tag}_{barcodes}_{AAorNT}-mutation-frequencies.csv',
        mutStats = '{tag}_mutation-stats.csv'
    output:
        'plots/{tag, [^\/_]*}_{barcodes, [^\/_]*}_{AAorNT, [^\/_]*}-mutation-frequencies.html'
    script:
        'utils/plot_mutations_frequencies.py'

rule hamming_distance:
    input:
        'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes.csv'
    output:
        HDmatrixCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA}-hamming-distance-matrix.csv',
        HDdistCSV = 'mutation_data/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_{NTorAA}-hamming-distance-distribution.csv'
    params:
        downsample = lambda wildcards: config.get('hamming_distance_distribution_downsample', False),
        refSeqs = lambda wildcards: config['runs'][wildcards.tag].get('reference', False)
    script:
        'utils/hamming_distance.py'

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
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False)
    script:
        'utils/plot_distribution.py'

rule plot_distribution_barcodeGroup:
    input:
        lambda wildcards: expand('mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_{NTorAA}-{distType}-distribution.csv',
                                tag = wildcards.tag,
                                barcodes = get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {})),
                                NTorAA = wildcards.NTorAA,
                                distType = wildcards.distType)
    output:
        plot = 'plots/{tag, [^\/_]*}_{NTorAA, [^\/_]*}-{distType, [^\/_]*}-distribution.html'
    params:
        labels = lambda wildcards: get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {})),
        title = lambda wildcards: wildcards.tag,
        legend_label = 'barcode group',
        background = lambda wildcards: config.get('background', False),
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False)
    script:
        'utils/plot_distribution.py'

rule plot_distribution_timepointGroup:
    input:
        lambda wildcards: expand('mutation_data/{tag_bc_tag_bc}_{NTorAA}-{distType}-distribution.csv',
                                tag_bc_tag_bc = [f'{tag}/{barcodes}/{tag}_{barcodes}' for tag,barcodes in get_demuxed_barcodes_timepoint( config['timepointsInfo'][wildcards.timepointsGroup]['tag_barcode_tp'].keys() )],
                                NTorAA = wildcards.NTorAA,
                                distType = wildcards.distType)
    output:
        plot = 'plots/timepoints/{timepointsGroup, [^\/_]*}_{NTorAA, [^\/_]*}-{distType, [^\/_]*}-distribution.html'
    params:
        labels = lambda wildcards: list(config['timepointsInfo'][wildcards.timepointsGroup]['tag_barcode_tp'].values()),
        title = lambda wildcards: wildcards.timepointsGroup,
        legend_label = lambda wildcards: config['timepointsInfo'][wildcards.timepointsGroup]['units'],
        background = lambda wildcards: config.get('background', False),
        raw = lambda wildcards: config.get('hamming_distance_distribution_raw', False)
    script:
        'utils/plot_distribution.py'

rule reduce_genotypes_dimensions:
    input:
        genotypes = '{dir}/{tag}_{barcodes}_genotypes.csv'
    output:
        reduced = '{dir}/{tag, [^\/_]*}_{barcodes, [^\/_]*}_genotypes-reduced-dimensions.csv'
    params:
        refSeqs = lambda wildcards: config['runs'][wildcards.tag].get('reference', False) if wildcards.tag in config['runs'] else config['timepointsInfo'][wildcards.tag].get('reference', False)
    script:
        'utils/dimension_reduction_genotypes.py'

rule plot_genotypes2D:
    input:
        genotypesReduced = 'mutation_data/{tag}/{barcodes}/{tag}_{barcodes}_genotypes-reduced-dimensions.csv'
    output:
        genotypes2Dplot = 'plots/{tag, [^\/_]*}/{barcodes, [^\/_]*}/{tag}_{barcodes}_genotypes2D.html'
    params:
        downsample = lambda wildcards: config.get('genotypes2D_plot_downsample', False),
        size_column = lambda wildcards: config.get('genotypes2D_plot_point_size_col', 'count'),
        size_range = lambda wildcards: config.get('genotypes2D_plot_point_size_range', '10, 30'),
        color_column = lambda wildcards: config.get('genotypes2D_plot_point_color_col', 'NT_substitutions_count')
    script:
        'utils/plot_genotypes_2d.py'

rule merge_timepoint_genotypes:
    input:
        genotypeCSVs = lambda wildcards: expand('mutation_data/{tag_barcodes}_genotypes.csv', tag_barcodes=[f'{tag}/{barcodes}/{tag}_{barcodes}' for tag,barcodes in get_demuxed_barcodes_timepoint( config['timepointsInfo'][wildcards.timepointsGroup]['tag_barcode_tp'].keys() )]),
    output:
        mergedGenotypes = 'mutation_data/timepoints/{timepointsGroup, [^\/_]*}_merged-timepoint-genotypes.csv'
    params:
        tpInfo = lambda wildcards: config['timepointsInfo'][wildcards.timepointsGroup]['tag_barcode_tp']
    run:
        import pandas as pd
        timepoints = list(params.tpInfo.values())   # timepoint values are an amount of time like X generations
        DFs = []
        for i, csv in enumerate(input.genotypeCSVs):
            df = pd.read_csv(csv)
            df['timepoint'] = timepoints[i]
            DFs.append(df)
        pd.concat(DFs).to_csv(output.mergedGenotypes)

rule reduce_genotypes_dimensions_timepoints:
    input:
        genotypes = 'mutation_data/timepoints/{timepointsGroup}_merged-timepoint-genotypes.csv'
    output:
        reduced = 'mutation_data/timepoints/{timepointsGroup, [^\/_]*}_merged-timepoint-genotypes-reduced-dimensions.csv'
    params:
        refSeqs = lambda wildcards: config['timepointsInfo'][wildcards.timepointsGroup]['reference']
    script:
        'utils/dimension_reduction_genotypes.py'

rule plot_genotypes2D_timepoints:
    input:
        genotypesReduced = 'mutation_data/timepoints/{timepointsGroup}_merged-timepoint-genotypes-reduced-dimensions.csv'
    output:
        genotypes2Dplot = 'plots/timepoints/{timepointsGroup, [^\/_]*}_genotypes2D.html'
    params:
        downsample = lambda wildcards: config.get('genotypes2D_plot_downsample', False),
        size_column = lambda wildcards: config.get('genotypes2D_plot_point_size_col', 'count'),
        size_range = lambda wildcards: config.get('genotypes2D_plot_point_size_range', '10, 30'),
        color_column = 'timepoint'
    script:
        'utils/plot_genotypes_2d.py'

def all_diversity_plots_input(wildcards):
    out = ( expand('plots/{tag}/{barcodes}/{tag}_{barcodes}_genotypes2D.html', tag=wildcards.tag, barcodes=get_demuxed_barcodes(wildcards.tag, config['runs'][wildcards.tag].get('barcodeGroups', {}))) )
    NTorAA = ['NT','AA'] if config['do_AA_mutation_analysis'][wildcards.tag] else ['NT']
    out.extend( expand('plots/{tag}_{NTorAA}-hamming-distance-distribution.html', tag=wildcards.tag, NTorAA=NTorAA) )
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
