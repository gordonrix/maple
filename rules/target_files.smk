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
    if any(config['do_UMI_analysis'][tag] for tag in config['runs']):
        out.append('sequences/UMI/UMI-extract-summary.csv')
        out.append('plots/UMI-group-distribution.html')
    if any(config['do_NT_mutation_analysis'][tag] for tag in config['runs']):
        out.append('mutation-stats.csv')
    if any(config['do_AA_mutation_analysis'][tag] for tag in config['runs']):
        if ('dms_view_chain' and 'dms_view_chain_numbering_difference') in config:
            out.append('dms-view-table.csv')
    if any(config['do_demux'][tag] for tag in config['runs']):
        out.append('demux-stats.csv')
        out.extend([f'plots/{tag}_demux.html' for tag in config['runs'] if config['do_demux'][tag]])
    for tag in config['runs']:
        if config['do_NT_mutation_analysis'][tag]:
            NTorAA = ['NT','AA'] if config['do_AA_mutation_analysis'][tag] else ['NT']
            plot_types = ['mutation-distribution-violin']
            if config.get('genotypes2D_plot_groups', False):
                plot_type.append('genotypes2D')
            out.extend(expand('plots/{tag}_{plot_type}.html', tag=tag, plot_type=plot_types))
            plot_type = ['mutation-distribution', 'mutation-frequencies', 'hamming-distance-distribution']
            out.extend(expand('plots/{tag}_{NTorAA}-{plot_type}.html', tag=tag, NTorAA=NTorAA, plot_type=plot_type))
            # out.extend(expand('plots/{tag}_mutation-spectra.html', tag=tag))
        if config['do_RCA_consensus'][tag]:
            out.append(f'plots/{tag}_RCA-distribution.html')
        if config['do_UMI_analysis'][tag]:
            if config['nanoplot'] == True:
                out.append(f"plots/nanoplot/{config['consensusCopyDict'][tag]}_alignment_preConsensus_NanoStats.txt")
        if config['nanoplot'] == True:
            out.append(f'plots/nanoplot/{tag}_fastq_NanoStats.txt')
            out.append(f'plots/nanoplot/{tag}_alignment_NanoStats.txt')
        if config['do_enrichment_analysis'].get(tag, False):
            out.append(f'plots/{tag}_enrichment-scores.html')            
        # out.append(f'plots/{tag}_pipeline-throughput.html')  # needs to be fixed to prevent use of temporary files that are computationally costly to recover

    # .done flag files are needed to separate the targets rule from rules that determine inputs from checkpoints (e.g. demux). A consequence of this is that
    #     simply deleting an input to one of these rules (what the user actually cares about) will not trigger the rule to be rerun. The user must delete
    #     the .done file as well (or just change a relevant parameter in the config) to trigger a rerun.
    if config.get('timepoints', {}):
        out.append('plots/.all_timepoints.done')

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

    if config.get('dashboard_input', False):
        db_input = dashboard_input(wildcards=None, config=config)
        if db_input:
            out.append(db_input['genotypes'])

    return out

rule targets:
    input:
        targets_input