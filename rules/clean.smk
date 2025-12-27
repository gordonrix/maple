# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : remove large files and move analysis files to timestamped folder
#                   so that analysis may be repeated
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------

# imports
import os, sys, glob
from datetime import datetime
# local rules
localrules: sequences_clean, alignment_clean, demux_clean, mutation_data_clean, logs_clean, clean

# Generate timestamp once and persist to file to ensure all rules use the same timestamp
# (Snakefile gets re-parsed multiple times, so datetime.now() would give different values)
timestamp_file = 'log/.clean_timestamp'
if not os.path.exists(timestamp_file):
    os.makedirs('log', exist_ok=True)
    timestamp_value = datetime.now().strftime("%Y%m%d-%H%M%S") + '-maple_outputs'
    with open(timestamp_file, 'w') as f:
        f.write(timestamp_value)
else:
    with open(timestamp_file, 'r') as f:
        timestamp_value = f.read().strip()

config['timestamp_dir'] = timestamp_value


# clean up sequences folder: remove combined .fastq.gz files (but not basecalled batches), move UMI stats files
rule sequences_clean:
    output:
        touch('log/.sequences_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp_dir']
    run:
        from utils.common import retrieve_fastqs
        import shutil

        # Create log directory if it doesn't exist
        os.makedirs('log', exist_ok=True)

        # Move UMI csv, fasta, and log files to timestamp directory
        if os.path.isdir('sequences/UMI'):
            umi_files = glob.glob('sequences/UMI/*.csv') + glob.glob('sequences/UMI/*.fasta*') + glob.glob('sequences/UMI/*.log')
            if umi_files:
                os.makedirs(os.path.join(params.timestamp_dir, 'sequences/UMI'), exist_ok=True)
                for f in umi_files:
                    shutil.move(f, os.path.join(params.timestamp_dir, 'sequences/UMI'))

        # Move RCA csv and fasta files to timestamp directory
        if os.path.isdir('sequences/RCA'):
            rca_files = glob.glob('sequences/RCA/*.csv') + glob.glob('sequences/RCA/*.fasta*')
            if rca_files:
                os.makedirs(os.path.join(params.timestamp_dir, 'sequences/RCA'), exist_ok=True)
                for f in rca_files:
                    shutil.move(f, os.path.join(params.timestamp_dir, 'sequences/RCA'))

        # Clean up paired-end failed merge files
        if os.path.isdir('sequences/paired'):
            for pattern in ['*_failed-merge_1.fastq.gz', '*_failed-merge_2.fastq.gz', '*_NGmerge.log']:
                for f in glob.glob(f'sequences/paired/{pattern}'):
                    os.remove(f)

        # Clean up sequences/{tag}.fastq.gz files after verifying originals still exist
        for tag in config['runs']:
            target_file = f'sequences/{tag}.fastq.gz'

            # Skip if target file doesn't exist
            if not os.path.exists(target_file):
                continue

            # Get original source files for this tag
            try:
                if 'fwdReads' in config['runs'][tag] or 'rvsReads' in config['runs'][tag]:
                    # Paired-end sequences
                    fwd_files = retrieve_fastqs(config['sequences_dir'], [config['runs'][tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][tag].get('fwdReads',''))
                    rvs_files = retrieve_fastqs(config['sequences_dir'], [config['runs'][tag].get('runname','')[0]], config['fastq_dir'], select=config['runs'][tag].get('rvsReads',''))
                    source_files = fwd_files + rvs_files
                elif 'runname' in config['runs'][tag]:
                    # Combined sequences
                    source_files = retrieve_fastqs(config['sequences_dir'], config['runs'][tag].get('runname',''), config['fastq_dir'])
                elif 'bs_project_ID' in config['runs'][tag]:
                    # BaseSpace sequences - skip verification as they may have been downloaded
                    print(f"[NOTICE] Skipping verification for BaseSpace sequences for tag {tag}")
                    continue
                else:
                    print(f"[WARNING] No source defined for tag {tag}, skipping")
                    continue

                # Verify source files exist
                missing_sources = [f for f in source_files if not os.path.exists(f)]
                if missing_sources:
                    print(f"[WARNING] Source files for {tag} not found, keeping sequences/{tag}.fastq.gz:")
                    for f in missing_sources:
                        print(f"  Missing: {f}")
                    continue

                # Get sizes of source and target files
                source_size = sum(os.path.getsize(f) for f in source_files)
                target_size = os.path.getsize(target_file)

                # Check if target is at least 90% of source size
                if target_size >= 0.9 * source_size:
                    os.remove(target_file)
                else:
                    size_pct = (target_size / source_size) * 100
                    print(f"[WARNING] sequences/{tag}.fastq.gz is only {size_pct:.1f}% of source size, keeping file")

            except Exception as e:
                print(f"[WARNING] Error verifying sources for {tag}: {e}")
                print(f"  Keeping sequences/{tag}.fastq.gz")

        # Clean up empty subdirectories in sequences folder, then sequences folder itself if empty
        for subdir in ['sequences/UMI', 'sequences/RCA', 'sequences/paired', 'sequences']:
            if os.path.isdir(subdir):
                contents = [f for f in os.listdir(subdir) if not f.startswith('.')]
                if not contents:
                    shutil.rmtree(subdir)
                    print(f"[NOTICE] Removed empty directory {subdir}")

        # Write output file
        with open(output[0], 'w') as f:
            f.write('done\n')

# clean up compute batches alignment
rule alignment_clean:
    output:
        touch('log/.alignment_clean.done')
    shell:
        """
        mkdir -p log
        if [ -d alignments ]; then
            rm -r alignments
        fi
        """

# clean up compute batches demux
rule demux_clean:
    input:
        'log/.mutation_data_clean.done'
    output:
        touch('log/.demux_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp_dir']
    shell:
        """
        mkdir -p log
        if [ -d demux ]; then
            rm -r demux
        fi
        """

rule mutation_data_clean:
    output:
        touch('log/.mutation_data_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp_dir'],
        keep = [directoryORfile for directoryORfile in os.listdir('.') if directoryORfile in ['plots', 'mutSpectra', 'mutation_data', 'enrichment', 'mutation-stats.csv', 'demux-stats.csv', 'dms-view-table.csv', 'maple']]
    shell:
        """
        mkdir -p log
        if [ ! -z "{params.keep}" ]; then
            mkdir -p {params.timestamp_dir}
            mv {params.keep} {params.timestamp_dir}/
        fi
        """

rule logs_clean:
    output:
        touch('log/.logs_clean.done')
    params:
        timestamp_dir = lambda wildcards: config['timestamp_dir'],
        metadata_dir = lambda wildcards: config['metadata']
    shell:
        """
        mkdir -p log
        log_files=(log/[!.]*)
        if [ -e "${{log_files[0]}}" ]; then
            mkdir -p {params.timestamp_dir}/mapleLogs
            mv "${{log_files[@]}}" {params.timestamp_dir}/mapleLogs
        fi
        snakemake_log_files=(.snakemake/log/[!.]*)
        if [ -e "${{snakemake_log_files[0]}}" ]; then
            mkdir -p {params.timestamp_dir}/snakemakeLogs
            mv "${{snakemake_log_files[@]}}" {params.timestamp_dir}/snakemakeLogs
        fi
        cp *.yaml {params.timestamp_dir}
        cp -r {params.metadata_dir} {params.timestamp_dir}
        """

# clean up everything
rule clean:
    input:
        rules.alignment_clean.output,
        rules.demux_clean.output,
        rules.mutation_data_clean.output,
        rules.logs_clean.output,
        rules.sequences_clean.output
    output:
        config['timestamp_dir'] + '.zip'
    params:
        timestamp_dir = lambda wildcards: config['timestamp_dir']
    shell:
        """
        rm {input}
        zip -r {params.timestamp_dir}.zip {params.timestamp_dir}
        rm -rf {params.timestamp_dir}
        rm -f log/.clean_timestamp
        """
