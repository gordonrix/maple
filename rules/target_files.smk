#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Defines rules whose sole purpose is to request outputs from other rules.
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

from utils.common import dashboard_input

def targets_input(wildcards):
    out = []

    if any(config['do_demux'][tag] for tag in config['runs']):
        out.append('demux-stats.csv')

    if any(config['do_UMI_analysis'][tag] for tag in config['runs']):
        out.append('sequences/UMI/UMI-extract-summary.csv')
        out.append('plots/UMI-group-distribution.html')
    if any(config['do_NT_mutation_analysis'][tag] for tag in config['runs']):
        out.append('mutation-stats.csv')
    if any(config['do_AA_mutation_analysis'][tag] for tag in config['runs']):
        if ('dms_view_chain' and 'dms_view_chain_numbering_difference') in config:
            out.append('dms-view-table.csv')

    for tag in config['runs']:

        if config['do_RCA_consensus'][tag]:
            out.append(f'plots/{tag}_RCA-distribution.html')
        if config['do_UMI_analysis'][tag]:
            if config['nanoplot'] == True:
                out.append(f"plots/nanoplot/{config['consensusCopyDict'][tag]}_alignment_preConsensus_NanoStats.txt")
        if config['nanoplot'] == True:
            out.append(f'plots/nanoplot/{tag}_fastq_NanoStats.txt')
            out.append(f'plots/nanoplot/{tag}_alignment_NanoStats.txt')

        # Enrichment outputs are now requested based on timepoint samples (see below)      

        #TODO: fix this      
        # out.append(f'plots/{tag}_pipeline-throughput.html')  # needs to be fixed to prevent use of temporary files that are computationally costly to recover

        if config['do_demux'][tag] and config['plot_demux']:
            out.append(f'plots/{tag}_demux.html')

        if config['do_NT_mutation_analysis'][tag]:
            NTorAA = ['NT','AA'] if config['do_AA_mutation_analysis'][tag] else ['NT']
            plot_types = []
            if config.get('mutations_violin_plot', False):
                plot_types.append('mutation-distribution-violin')
            if config.get('genotypes2D_plot_groups', False):
                plot_type.append('genotypes2D')
            plot_types = [PT for PT in plot_types if config.get(f'plot_{PT}', False)]
            out.extend(expand('plots/{tag}_{plot_type}.html', tag=tag, plot_type=plot_types))
            NTAA_plot_types = ['mutation-distribution', 'mutations-aggregated']
            if config.get('hamming_distance_distribution', False):
                NTAA_plot_types.append('hamming-distance-distribution')
            NTAA_plot_types = [PT for PT in NTAA_plot_types if config.get(f'plot_{PT}', False)]
            out.extend(expand('plots/{tag}_{NTorAA}-{plot_type}.html', tag=tag, NTorAA=NTorAA, plot_type=NTAA_plot_types))

    # .done flag files are needed to separate the targets rule from rules that determine inputs from checkpoints (e.g. demux). A consequence of this is that
    #     simply deleting an input to one of these rules (what the user actually cares about) will not trigger the rule to be rerun. The user must delete
    #     the .done file as well (or just change a relevant parameter in the config) to trigger a rerun.

    if config.get('timepointsInfo', {}):
        out.append('plots/.all_timepoints.done')

        # Request enrichment score outputs if enrichment is enabled
        if config.get('do_enrichment', False):
            enrichment_type = config.get('enrichment_type', 'genotype')
            timepoint_samples = list(config['timepointsInfo'].keys())

            # Request individual replicate scores for all samples
            out.extend(expand('enrichment/{sample}_{mode}-enrichment-scores.csv',
                             sample=timepoint_samples,
                             mode=enrichment_type))

            # Request mean scores for all groups (including single-replicate groups)
            # For single-replicate groups, the "mean" is just that replicate's scores
            groups = list(set([
                config['timepointsInfo'][sample]['group']
                for sample in timepoint_samples
            ]))

            out.extend(expand('enrichment/{group}_{mode}-enrichment-scores-mean.csv',
                             group=groups,
                             mode=enrichment_type))

            # Request dimension-reduced enrichment files for dashboard
            out.extend(expand('mutation_data/timepoints/{group}_merged-timepoint_genotypes-reduced-dimensions-{mode}-enrichment-mean.csv',
                             group=groups,
                             mode=enrichment_type))

    if config.get('genotypes2D_plot_all', False):
        out.append('plots/.all_genotypes_2D.done')

    if config.get('diversity_plot_subset', False):
        for tag_bc in config['diversity_plot_subset'].split(','):
            divPlotFilePrefixes = []
            plotType = ['genotypes2D.html', 'NT-hamming-distance-distribution.html']
            if config['do_AA_mutation_analysis'][tag]:
                plotType.append('AA-hamming-distance-distribution.html')
            tag, bc = tag_bc.split('_')
            divPlotFilePrefixes.append(f'{tag}/{bc}/{tag_bc}')
            out.extend( expand('plots/{tag_barcodes}_{plotType}', tag_barcodes=divPlotFilePrefixes, plotType=plotType) )

    out_string = "\n" + "\n".join(out) + "\n" # for debugging

    if config.get('dashboard_input', False):
        db_input = dashboard_input(wildcards=None, config=config)
        if db_input:
            out.append(db_input['datasets_csv'])
    return out

rule targets:
    input:
        targets_input