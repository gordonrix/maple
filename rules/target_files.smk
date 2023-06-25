#
#  DESCRIPTION   : Supplementary snakefile in the Maple snakemake pipeline.
#                   Defines rules whose sole purpose is to request outputs from other rules.
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

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
            plot_type = ['genotypes2D', 'mutation-distribution-violin']
            out.extend(expand('plots/{tag}_{plot_type}.html', tag=tag, plot_type=plot_type))
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
        if config.get('timepoints', {}).get(tag, False) and config['do_NT_mutation_analysis'][tag] == True:
            out.extend(expand('plots/{tag}_mutation-rates-mut-grouped.html', tag=config['timepoints']))                 # mutation rates includes data for all rows in the timepoints file
            out.extend(expand('plots/timepoints/{timepointSample}_{plot_type}.html', timepointSample=config['timepointsInfo'], plot_type=['genotypes2D', 'mutation-distribution-violin']))    # each of these outputs includes data for a single row in the timepoints file
            timepointSample_NTorAA = []
            out.extend(expand('plots/timepoints/{timepointSample_NTorAA}-{distType}-distribution.html', timepointSample_NTorAA=timepointSample_NTorAA, distType=['mutation', 'hamming-distance']))
            for timepointSample in config['timepointsInfo']:
                tag = config['timepointsInfo'][timepointSample]['tag']
                timepointSample_NTorAA.extend(expand('{TS}_{NTorAA}', TS=timepointSample, NTorAA= ['NT','AA'] if config['do_AA_mutation_analysis'][ tag ] else ['NT']))
        if config['do_enrichment_analysis'].get(tag, False):
            out.append(f'plots/{tag}_enrichment-scores.html')
            
        # out.append(f'plots/{tag}_pipeline-throughput.html')  # needs to be fixed to prevent use of temporary files that are computationally costly to recover
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

    return out

rule targets:
    input:
        targets_input