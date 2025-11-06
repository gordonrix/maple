#
#  DESCRIPTION   : Main snakefile for the Maple pipeline. Performs quality control and imports modules that
#                   perform analysis steps
#
#  RESTRICTIONS  : none
#
#  AUTHOR(S)     : Gordon Rix
#

# imports
import os, sys, collections, glob
import itertools
import yaml, subprocess
import snakemake.common
from datetime import datetime
from snakemake.utils import min_version
from Bio import SeqIO
import json
import yaml
import pandas as pd
import numpy as np
import logging
from rules.utils.common import load_csv_as_dict, validate_bc, str_to_bool

start_time = datetime.now()


# snakemake config
min_version("5.5.2")
configfile: "config.yaml"


# get pipeline version # update for maple
def get_tag():
    try:
        cmd = 'git describe --tags'
        version = subprocess.check_output(cmd.split(), cwd=os.path.dirname(workflow.snakefile)).decode().strip()
    except subprocess.CalledProcessError:
        print('[WARNING] Unable to get version from git tags.', file=sys.stderr)
        version = '-'
    try:
        cmd = 'git rev-parse --abbrev-ref HEAD'
        branch = subprocess.check_output(cmd.split(), cwd=os.path.dirname(workflow.snakefile)).decode().strip()
    except subprocess.CalledProcessError:
        print('[WARNING] Unable to get branch from git. Pulling development.', file=sys.stderr)
        branch = 'development'
    if '-' in version:
        if branch == 'master':
            return 'latest', version
        else:
            return branch, version
    else:   # clean tag checkout
        return version, version


maple_tag, maple_git_tag = get_tag()
config['version'] = {'tag': maple_tag, 'full-tag': maple_git_tag}


# scan working directory
def get_dir_files(base_dir):
    return set({os.path.join(path, name) for path, subdirs, files in os.walk(base_dir) for name in files if not '/.' in path})

start_files = get_dir_files(workflow.workdir_init)

# add memory scaling for minimap2
config['memory'] = {}
config['memory']['minimap2'] = [8000,500]

# list of required options. Will be added to if certain tools are used
required = ['sequences_dir', 'fastq_dir', 'metadata', 'threads_alignment', 'threads_samtools', 'threads_demux', 'NGmerge_flags', 'nanoplot', 'nanoplot_flags', 'alignment_samtools_flags', 'alignment_minimap2_flags', 'mutation_analysis_quality_score_minimum', 'sequence_length_threshold', 'highest_abundance_genotypes', 'mutations_frequencies_raw', 'unique_genotypes_count_threshold', 'runs']
if ('minknowDir' in config) and ('sequences_dir' not in config):
    print('[WARNING] "sequences_dir" is a new required config variable that replaces "minknowDir", please replaced "minknowDir" with "sequences_dir" in your config file"')

def str_to_number(value):
    """
    Convert a string to a number (int or float) if possible, raise error if not.
    """
    try:
        # Try converting to integer first
        number = int(value)
        return number
    except ValueError:
        raise ValueError(f"Invalid number string: {value}")
    return number


def pretty_print(d, indent=2):
    """
    Recursively prints a nested dictionary with indentation.
    
    :param d: The dictionary to print
    :param indent: The current level of indentation (used internally)
    """
    for key, value in d.items():
        print('  ' * indent + str(key) + ":")
        if isinstance(value, dict):
            pretty_print(value, indent + 1)
        else:
            print('  ' * (indent + 1) + str(value) + ' ' + type(value).__name__)

def validate_csv_path(csv_path, metadata_dir, config_key, tag):
    """
    Validate that a config value is a string path to a CSV file.
    Returns tuple: (full_path, error_message or None)
    """
    if not isinstance(csv_path, str):
        return None, f"[ERROR] {config_key} for tag `{tag}` must be a string path to a CSV file, not {type(csv_path).__name__}.\n"
    
    # Add metadata directory if not already present
    if not csv_path.startswith(metadata_dir):
        full_path = os.path.join(metadata_dir, csv_path)
    else:
        full_path = csv_path
    
    # Check if file exists and is a CSV
    if not os.path.isfile(full_path):
        return None, f"[ERROR] {config_key} file `{full_path}` for tag `{tag}` does not exist. {csv_path}\n"
    
    if not full_path.endswith('.csv'):
        return None, f"[ERROR] {config_key} file `{full_path}` for tag `{tag}` does not have .csv extension.\n"
    
    return full_path, None


def validate_barcode_groups_csv(bc_groups_csv, barcode_info_csv, tag, group_type='partition'):
    """
    Validate barcode groups CSV file (either partition or label groups).
    Returns a list of warnings.
    """
    warnings = []

    try:
        groups_df = pd.read_csv(bc_groups_csv)
    except Exception as e:
        warnings.append(f"[WARNING] Failed to read {group_type} barcode groups CSV file `{bc_groups_csv}` for tag `{tag}`: {e}\n")
        return warnings

    # Filter by tag if tag column exists
    if 'tag' in groups_df.columns:
        groups_df = groups_df[groups_df['tag'] == tag]
    if 'barcode_group' not in groups_df.columns:
        warnings.append(f"[WARNING] Required column `barcode_group` not found in {group_type} groups CSV `{bc_groups_csv}` for tag `{tag}`.\n")
        return warnings

    if len(groups_df) == 0:
        return warnings  # Empty groups are OK

    # Load barcode info to check against
    try:
        barcode_info_df = pd.read_csv(barcode_info_csv)
        if 'tag' in barcode_info_df.columns:
            barcode_info_df = barcode_info_df[barcode_info_df['tag'] == tag]

        # Create barcode info lookup
        barcode_info = {}
        for _, row in barcode_info_df.iterrows():
            barcode_info[row['barcode_name']] = {
                'label_only': row.get('label_only', False),
                'fasta': row['fasta']
            }
    except Exception as e:
        warnings.append(f"[WARNING] Could not read barcode info CSV for validation of {group_type} groups: {e}\n")
        return warnings

    barcode_types = list(groups_df.columns.difference(['barcode_group', 'tag']))

    # Check each group
    for _, row in groups_df.iterrows():
        group_name = row['barcode_group']

        if '_' in group_name:
            warnings.append(f"[WARNING] {group_type.capitalize()} barcode group `{group_name}` for run tag `{tag}` contains underscore(s), which will disrupt the pipeline.\n")

        # Check each barcode type in the group
        for col in barcode_types:

            if pd.notna(row[col]):
                barcode_type = col
                barcode_name = row[col]

                # Check if barcode type exists
                if barcode_type not in barcode_info:
                    warnings.append(f"[WARNING] Barcode type `{barcode_type}` in {group_type} group `{group_name}` for run tag `{tag}` is not defined in barcode info.\n")
                    continue

                # Check label_only status
                is_label_only = barcode_info[barcode_type].get('label_only', False)
                is_label_only = False if np.isnan(is_label_only) else is_label_only

                if group_type == 'partition' and is_label_only:
                    warnings.append(f"[WARNING] Barcode type `{barcode_type}` is marked as label_only but is used in partition group `{group_name}` for tag `{tag}`. Demultiplexing will fail.\n")
                elif group_type == 'label' and not is_label_only:
                    warnings.append(f"[WARNING] Barcode type `{barcode_type}` is not marked as label_only but is used in label group `{group_name}` for tag `{tag}`. Demultiplexing will fail.\n")

    return warnings


def validate_generated_barcode_consistency(config):
    """
    Validate that barcode types with generate=True have consistent settings
    across all tags that share the same output fasta file.
    Returns a list of error messages.
    """
    errors = []

    if 'barcode_fasta_to_tags' not in config:
        return errors

    fasta_to_barcode_settings = {}

    for fasta_path, tags_list in config['barcode_fasta_to_tags'].items():
        for tag in tags_list:
            barcode_info = config['runs'][tag]['barcode_info']

            for barcode_name, info in barcode_info.items():
                # Only check barcode types with generate=True
                if not info.get('generate', False):
                    continue

                # Check if this barcode type outputs to this fasta file
                bc_fasta_path = info['fasta']
                if not bc_fasta_path.startswith(config['metadata']):
                    bc_fasta_path = os.path.join(config['metadata'], bc_fasta_path)

                if bc_fasta_path != fasta_path:
                    continue

                # Store or validate settings for this fasta file
                if fasta_path not in fasta_to_barcode_settings:
                    fasta_to_barcode_settings[fasta_path] = {
                        'barcode_name': barcode_name,
                        'context': info['context'],
                        'hamming_distance': info.get('hamming_distance', 0),
                        'reverse_complement': info.get('reverse_complement', False),
                        'first_tag': tag
                    }
                else:
                    # Validate consistency
                    stored = fasta_to_barcode_settings[fasta_path]
                    mismatches = []

                    if stored['context'] != info['context']:
                        mismatches.append(f"context: '{stored['context']}' vs '{info['context']}'")
                    if stored['hamming_distance'] != info.get('hamming_distance', 0):
                        mismatches.append(f"hamming_distance: {stored['hamming_distance']} vs {info.get('hamming_distance', 0)}")
                    if stored['reverse_complement'] != info.get('reverse_complement', False):
                        mismatches.append(f"reverse_complement: {stored['reverse_complement']} vs {info.get('reverse_complement', False)}")

                    if mismatches:
                        error_msg = (
                            f"[ERROR] Inconsistent barcode settings for generated fasta file '{fasta_path}':\n"
                            f"  Tag '{stored['first_tag']}' (barcode '{stored['barcode_name']}') vs "
                            f"tag '{tag}' (barcode '{barcode_name}'):\n"
                            f"  {', '.join(mismatches)}\n"
                            f"  All tags sharing a barcode fasta file with generate=True must have identical settings "
                            f"for context, hamming_distance, and reverse_complement.\n"
                        )
                        errors.append(error_msg)

    return errors

# keep a list of errors, warnings, and notices
errors = []
warnings = []
notices = []
# notify user if any notices were not sent to stdout
notice_of_notices = False
runs_to_import = []


# import and validate any csv metadata files

if type(config.get('runs', False)) is str:
    full_csv_path = os.path.join(config['metadata'], config['runs'])
    df = pd.read_csv(full_csv_path, index_col=False)
    config['runs'] = load_csv_as_dict(full_csv_path, required=['reference'], lists=['runname', 'UMI_contexts'])

for tag in config['runs']:
    # Check for deprecated barcodeInfo dictionary
    if 'barcodeInfo' in config['runs'][tag]:
        if isinstance(config['runs'][tag]['barcodeInfo'], dict):
            errors.append(f"[ERROR] Dictionary input for barcodeInfo for tag `{tag}` is no longer supported. Please provide a CSV file path with the parameter `barcode_info_csv` instead.\n")
    
    # Check for deprecated barcodeGroups dictionary
    if 'barcodeGroups' in config['runs'][tag]:
        if isinstance(config['runs'][tag]['barcodeGroups'], dict):
            errors.append(f"[ERROR] Dictionary input for barcodeGroups for tag `{tag}` is no longer supported. Please provide a CSV file path with the parameter `partition_barcode_groups_csv` instead.\n")

    # Handle barcode_info_csv (validate and convert to full path)
    barcode_info_csv = config['runs'][tag].get('barcode_info_csv', None)
    if isinstance(barcode_info_csv, str):
        csv_path, error = validate_csv_path(
            config['runs'][tag]['barcode_info_csv'], 
            config['metadata'], 
            'barcode_info_csv',
            tag
        )
        if error:
            errors.append(error)
        else:
            config['runs'][tag]['barcode_info_csv'] = csv_path

        # add barcode info as dictionary to config
        config['runs'][tag]['barcode_info'] = load_csv_as_dict(
            csv_path,
            required=['barcode_name', 'context', 'fasta', 'reverse_complement'],
            tag=tag,
            defaults={'reverse_complement': False, 'label_only': False, 'generate': False, 'hamming_distance': 0}
        )

        # Build barcode fasta to tags mapping (for generating shared barcode references)
        if 'barcode_fasta_to_tags' not in config:
            config['barcode_fasta_to_tags'] = {}

        barcode_info = config['runs'][tag]['barcode_info']
        for barcode_name, info in barcode_info.items():
            if info.get('generate'):
                fasta_path = info['fasta']
                if not fasta_path.startswith(config['metadata']):
                    fasta_path = os.path.join(config['metadata'], fasta_path)
                if fasta_path not in config['barcode_fasta_to_tags']:
                    config['barcode_fasta_to_tags'][fasta_path] = []
                if tag not in config['barcode_fasta_to_tags'][fasta_path]:
                    config['barcode_fasta_to_tags'][fasta_path].append(tag)

        # Clean up barcode_info: remove empty entries, ensure inputs adhere to expected values
        for barcode_name in barcode_info:

            config['runs'][tag]['barcode_info'], bc_errors = validate_bc(
                config['runs'][tag]['barcode_info'], 
                barcode_name,
                metadata_dir=config['metadata'],
                errors_list=errors  # Pass the global errors list
            )

    # Handle partition_barcode_groups_csv and label_barcode_groups_csv (validate and convert to full path)

    for group_type in ['partition', 'label']:
        key = f'{group_type}_barcode_groups_csv'
        csv = config['runs'][tag].get(key, None)
        if isinstance(csv, str):
            csv_full_path, error = validate_csv_path(
                csv,
                config['metadata'],
                key,
                tag
            )
            if error:
                errors.append(error)
            else:
                config['runs'][tag][key] = csv_full_path

            # Validate groups if barcode_info_csv exists
            if 'barcode_info_csv' in config['runs'][tag]:
                val_warnings = validate_barcode_groups_csv(
                    csv_full_path,
                    config['runs'][tag]['barcode_info_csv'],
                    tag,
                    group_type=group_type
                )
                for w in val_warnings:
                    print(w, file=sys.stderr)
            else:
                warnings.append(f"[WARNING] {key} provided but no barcode_info_csv provided for tag `{tag}`. Demultiplexing will not be run.\n")

# Validate consistency of barcode settings across tags that share generated barcode fasta files
barcode_consistency_errors = validate_generated_barcode_consistency(config)
errors.extend(barcode_consistency_errors)

# check for sequences
for tag in config['runs']:
    sequences = os.path.join('sequences', tag+'.fastq.gz')

    if not os.path.exists(sequences):

        if any([x in config['runs'][tag] for x in ['fwdReads', 'rvsReads']]):
            if 'fwdReads' not in config['runs'][tag] or 'rvsReads' not in config['runs'][tag]:
                print(f"[WARNING] Both forward and reverse reads files not provided for {tag} with keyword `fwdReads` and `rvsReads`.\n", file=sys.stderr)
            fwd = os.path.join('sequences', 'paired', config['runs'][tag]['fwdReads'])
            rvs = os.path.join('sequences', 'paired', config['runs'][tag]['rvsReads'])
            if not all((os.path.exists(fwd), os.path.exists(rvs))):
                print(f"[WARNING] One or both forward/reverse reads files provided for {tag}, ({fwd}, {rvs}) do not exist.\n", file=sys.stderr)

        elif 'runname' in config['runs'][tag]:
            for runname in config['runs'][tag]['runname']:
                batch = os.path.join('sequences', 'batches', runname)
                if not os.path.exists(batch):
                    runs_to_import.append(runname)

        read_source_key_count = sum(1 for key in ['fwdReads', 'runname', 'bs_project_ID'] if key in config['runs'][tag])
        if read_source_key_count > 1:
            errors.append(f"[ERROR] Multiple sources for reads provided for tag `{tag}`. Please provide only one of (`fwdReads`/`rvsReads`), `runname`, or `bs_project_ID`.\n")

# Check sequences directory
if runs_to_import != []:
    for runname in runs_to_import:
        if not os.path.isdir(config['sequences_dir']):
            print(f"[WARNING] May need to import runname `{runname}`, but the provided sequences directory, `{config['sequences_dir']}`, does not exist.\n", file=sys.stderr)
        else:
            # Check if runname is a fastq file
            is_fastq = any(runname.endswith(ext) for ext in ['.fastq', '.fastq.gz', '.fq.gz'])
            
            if is_fastq:
                # Search for the specific file
                import glob
                pattern = os.path.join(config['sequences_dir'], '**', runname)
                matches = glob.glob(pattern, recursive=True)
                
                if not matches:
                    print(f"[WARNING] May need to import fastq file `{runname}`, but this could not be located in any directory tree under `{config['sequences_dir']}`.\n", file=sys.stderr)
                elif len(matches) > 1:
                    print(f"[WARNING] More than one ({matches}) of fastq file `{runname}` found under `{config['sequences_dir']}`\n", file=sys.stderr)
            else:
                # Original folder check
                if all([runname not in dirs for _, dirs, _ in os.walk(config['sequences_dir'].rstrip('/'))]):
                    print(f"[WARNING] May need to import runname `{runname}`, but this could not be located in any directory tree under `{config['sequences_dir']}`.\n", file=sys.stderr)


# check reference sequences
refSeqFastaFiles = []   # list files for all tags then check so files not being checked multiple times
config['do_NT_mutation_analysis'] = {}     # dictionaries to determine if NT/AA analysis should be performed based on how many sequences are present in the ref fasta file
config['do_AA_mutation_analysis'] = {}
for tag in config['runs']:
    if 'reference' not in config['runs'][tag]:
        errors.append(f"[ERROR] No reference file provided for tag `{tag}")
    refPath = config['runs'][tag]['reference']
    refName = os.path.basename(refPath)
    refFullPath = os.path.join(config['metadata'], refPath)
    ref_parent_dir = os.path.dirname(refFullPath)
    if not (refName.endswith('fasta') or refName.endswith('.fa')):
        print(f'[WARNING] Reference .fasta file for {tag} does not end with `.fasta` or `.fa` (given path: {refFullPath}).', file=sys.stderr)
    config['runs'][tag]['reference'] = refFullPath
    if not os.path.isfile(refFullPath):
        print(f'[ERROR] Reference .fasta file for {tag} (given path: {refFullPath}) not found.', file=sys.stderr)
    refFastaPrefix = refName.split('.f')[0]
    alnRefFullPath = os.path.join(config['metadata'], '.' + refFastaPrefix + '_aln.fasta')
    config['runs'][tag]['reference_aln'] = alnRefFullPath
    if (refFullPath, alnRefFullPath) not in refSeqFastaFiles:
        refSeqFastaFiles.append((refFullPath, alnRefFullPath))
    referenceSeqs = list(SeqIO.parse(refFullPath, 'fasta'))

    config['do_NT_mutation_analysis'][tag] = False
    config['do_AA_mutation_analysis'][tag] = False
    if len(referenceSeqs) >= 2:
        config['do_NT_mutation_analysis'][tag] = True
    if len(referenceSeqs) == 3:
        config['do_AA_mutation_analysis'][tag] = True
    if ('AA_muts_of_interest' in config['runs'][tag]) and not config['do_AA_mutation_analysis'][tag]:
        print(f'[WARNING] AA_muts_of_interest provided for run tag `{tag}`, but no protein seq provided for this tag. AA muts of interest will not be evaluated for this tag.', file=sys.stderr)

for refFasta, alnFasta in refSeqFastaFiles:
    referenceSeqs = list(SeqIO.parse(refFasta, 'fasta'))

    if len(referenceSeqs) not in [1,2,3]:
        errors.append(f"[ERROR] Reference sequence file {refFasta} Does not contain 1, 2, or 3 sequences. Ensure file is fasta formatted and does not contain erroneous sequences.")
    alignmentSeq, nucleotideSeq, proteinSeq = False, False, False
    alignmentSeq = referenceSeqs[0]
    if len(referenceSeqs) >= 2:
        nucleotideSeq = referenceSeqs[1]
    if len(referenceSeqs) == 3:
        proteinSeq = referenceSeqs[2]

    if nucleotideSeq:
        if alignmentSeq.seq.upper().find(nucleotideSeq.seq.upper()) == -1:
            if alignmentSeq.seq.upper().find(nucleotideSeq.seq.reverse_complement().upper()) == -1:
                errors.append(f"[ERROR] Nucleotide (second) sequence, `{nucleotideSeq.id}`, nor its reverse complement, is not a subsequence of alignment (first) sequence, `{alignmentSeq.id}`, in reference file `{refFasta}`.\n")

    if proteinSeq:
        if nucleotideSeq.seq.upper().find(proteinSeq.seq.upper()) == -1:
            errors.append(f"[ERROR] Protein (third) sequence, `{proteinSeq.id}`, is not a subsequence of nucleotide (second) sequence, {nucleotideSeq.id}, in reference file `{refFasta}`.\n")
        if len(proteinSeq.seq)%3 != 0:
            errors.append(f"[ERROR] Length of protein reference sequence `{proteinSeq.id}` of reference file `{refFasta}` is not a multiple of 3, and therefore cannot be used as ORF\n")
        for i, nt in enumerate(str(proteinSeq.seq).upper()):
            if nt not in list("ATGCN"):
                errors.append(f"[ERROR] Character {nt} at position {i} in reference sequence `{proteinSeq.id}` of reference file `{refFasta}` is not a canonical nucleotide or 'N'\n")

    # auto generate file used for alignment so that cropping / extending other sequences(es) in refFasta doesn't command a re-run of time consuming steps like alignment and UMI consensus generation
    if os.path.isfile(alnFasta):
        try:
            refFirstRecord = next(SeqIO.parse(refFasta, 'fasta'))
            alnFirstRecord = next(SeqIO.parse(alnFasta, 'fasta'))
            refFirstRecord.seq = refFirstRecord.seq.upper()
            alnFirstRecord.seq = alnFirstRecord.seq.upper()
            # make new file if aln record not the same as first record from ref
            if (refFirstRecord.seq != alnFirstRecord.seq) or (refFirstRecord.id != alnFirstRecord.id):
                os.remove(alnFasta)
        except StopIteration:
            os.remove(alnFasta)
    if not os.path.isfile(alnFasta):
        print(f'Alignment reference .fasta file not found or is different from original reference .fasta file. Generating {alnFasta} from {refFasta}.\n', file=sys.stderr)
        with open(alnFasta, 'w') as fastaOut:
            first_record = next(SeqIO.parse(refFasta, 'fasta'))
            fastaOut.write(f'>{first_record.id}\n{first_record.seq.upper()}\n')

# RCA/UMI consensus checks and consensus copy dict
config['do_RCA_consensus'] = {}
config['do_UMI_analysis'] = {}
consensusCopyDict = {}      # dictionary that is used to reuse consensus sequences from a different tag if they are generated using the same files. Keys are tags, and values are tags whose consensus sequences will be used for downstream files for the key tag
consensusRecipeDict = {}     # dictionary that keeps track of the runnames, reference sequence, UMI contexts, and splint sequence used to generate a consensus for each tag. If any two tags share all of these, they will use the same consensus sequence output to conserve computation and storage
for tag in config['runs']:

    if 'splint' in config['runs'][tag]:
        config['do_RCA_consensus'][tag] = True
        splintSeq = config['runs'][tag]['splint']
        splintFasta = os.path.join(config['metadata'], f".{tag}_splint.fasta")
        makeSplintFasta = True
        if os.path.isfile(splintFasta):
            first_record = next(SeqIO.parse(splintFasta, 'fasta'))
            if (first_record.id == 'splint') and (str(first_record.seq).upper() == splintSeq.upper()):
                makeSplintFasta = False
            else:
                print(f"[NOTICE] Splint sequence provided in config file for run tag `{tag}` has changed. Updating splint fasta file `{splintFasta}`.\n", file=sys.stderr)
        if makeSplintFasta:
            with open(splintFasta, 'w') as out:
                out.write(f'>splint\n{splintSeq}')
    else:
        config['do_RCA_consensus'][tag] = False

    refFasta = config['runs'][tag]['reference']
    alignmentSeq = list(SeqIO.parse(refFasta, 'fasta'))[0]
    if 'UMI_contexts' in config['runs'][tag]:
        config['do_UMI_analysis'][tag] = True
        config['runs'][tag]['UMI_contexts'] = [context.upper() for context in config['runs'][tag]['UMI_contexts']]
        if len(set(config['runs'][tag]['UMI_contexts'])) != len(config['runs'][tag]['UMI_contexts']):
            errors.append(f"[ERROR] Duplicate UMI contexts provided for tag `{tag}`. UMI consensus generation will fail.\n")
        for i, context in enumerate(config['runs'][tag]['UMI_contexts']):
            occurences = str(alignmentSeq.seq).upper().count(context.upper())
            if occurences == 0:
                errors.append(f"[ERROR] UMI context {i+1} for tag `{tag}`, `{context}`, not found in reference `{alignmentSeq.id}` in fasta `{refFasta}`. UMI consensus generation will fail.\n")
            elif occurences > 1:
                errors.append(f"[ERROR] UMI context {i+1} for tag `{tag}`, `{context}`, present more than once in reference `{alignmentSeq.id}` in fasta `{refFasta}`. UMI consensus generation will fail.\n")
    else:
        config['do_UMI_analysis'][tag] = False

    if ( config['do_UMI_analysis'][tag] ) or ( config['do_RCA_consensus'][tag] ):
        
        if 'runname' in config['runs'][tag]:
            rawdata = config['runs'][tag]['runname']
            rawdata.sort()
            rawdata = tuple(rawdata)
        elif 'fwdReads' in config['runs'][tag]:
            rawdata = (config['runs'][tag]['fwdReads'], config['runs'][tag]['rvsReads'])
        else:
            rawdata = tag   # if using manually input data, dont use for multiple tags

        UMIcontexts = config['runs'][tag].get('UMI_contexts', [])
        UMIcontexts.sort()
        UMIcontexts = tuple(UMIcontexts)

        consensusRecipe = ( rawdata, UMIcontexts, config['runs'][tag].get('splint', ''), str(alignmentSeq.seq).upper() )

        if consensusRecipe in consensusRecipeDict:
            print(f"[NOTICE] Raw data / alignment sequence / UMI context / splint combination used more than once. Using the consensus .fasta file of tag `{consensusRecipeDict[consensusRecipe]}` for tag `{tag}` to reduce computation time and storage requirements.\n", file=sys.stderr)
        else:
            consensusRecipeDict[consensusRecipe] = tag
        consensusCopyDict[tag] = consensusRecipeDict[consensusRecipe]

if any([config['do_UMI_analysis'][tag] for tag in config['runs']]):
    required.extend(['medaka_model', 'medaka_flags', 'UMI_medaka_batches', 'threads_medaka', 'UMI_mismatches', 'UMI_consensus_minimum', 'UMI_consensus_maximum'])
if any([config['do_RCA_consensus'][tag] for tag in config['runs']]):
    required.extend(['peak_finder_settings', 'RCA_batch_size', 'RCA_consensus_minimum', 'RCA_consensus_maximum'])

config['consensusCopyDict'] = consensusCopyDict

# check that any global config options that use dictionaries contain all tags as keys
def validate_params(potential_param_dicts):
    param_dicts = [param for param in potential_param_dicts if type(config[param])==dict]
    for tag in config['runs']:
        for param in param_dicts:
            if tag not in config[param].keys():
                print(f"[WARNING] Dict used for parameter {param}, but run tag {tag} not found as a key in this dict.\n")
    return

potential_param_dicts = ['alignment_minimap2_flags']
validate_params(potential_param_dicts)

# Demultiplexing checks
config['do_demux'] = {}
for tag in config['runs']:
    if 'barcode_info_csv' not in config['runs'][tag]:
        config['do_demux'][tag] = False
        continue
    else:
        config['do_demux'][tag] = True
        barcode_info = config['runs'][tag]['barcode_info']
    
    # Check partition groups
    if 'partition_barcode_groups_csv' in config['runs'][tag]:
        try:
            groups_df = pd.read_csv(config['runs'][tag]['partition_barcode_groups_csv'])
            if 'tag' in groups_df.columns:
                groups_df = groups_df[groups_df['tag'] == tag]
            if len(groups_df) == 0:
                print(f"[NOTICE] No partition barcode groups provided for run tag `{tag}`. Will use barcode names to name demultiplexed files.\n", file=sys.stderr)
        except Exception:
            pass
    else:
        notice = f"[NOTICE] `barcode_info_csv` supplied but `partition_barcode_groups_csv` not supplied for run tag `{tag}`. Will use barcode names to name demultiplexed files.\n"
        notices.append(notice)
        notice_of_notices = True
    
    # Check label groups
    if 'label_barcode_groups_csv' in config['runs'][tag]:
        try:
            label_df = pd.read_csv(config['runs'][tag]['label_barcode_groups_csv'])
            if 'tag' in label_df.columns:
                label_df = label_df[label_df['tag'] == tag]
            if len(label_df) == 0:
                print(f"[NOTICE] No label barcode groups provided for run tag `{tag}`. Will use barcode names to label sequences.\n", file=sys.stderr)
        except Exception:
            pass
    
    # Get reference sequence
    refFasta = config['runs'][tag]['reference']
    alignmentSeq = list(SeqIO.parse(refFasta, 'fasta'))[0]
    contexts = {'partition': [], 'label': []}  # contexts for partition and label barcodes
    
    # Validate each barcode type
    for barcodeType in barcode_info:
        barcode_data = barcode_info[barcodeType]
        
        # Check context
        c_given = barcode_data.get('context', False)
        context = c_given.upper()
        label_only = barcode_data.get('label_only', False)
        
        if context:
            if label_only:
                contexts['label'].append(context)
            else:
                contexts['partition'].append(context)
        
        # Check context in reference
        occurrences = str(alignmentSeq.seq).upper().count(context)
        if occurrences == 0:
            print(f"[WARNING] Barcode type `{barcodeType}` context `{c_given}` not found in reference `{alignmentSeq.id}` in fasta `{refFasta}`\n", file=sys.stderr)
        elif occurrences > 1:
            print(f"[WARNING] Barcode type `{barcodeType}` context `{c_given}` found more than once in reference `{alignmentSeq.id}` in fasta `{refFasta}`\n", file=sys.stderr)
        
        # Check barcode fasta
        bcFasta = barcode_data.get('fasta', '')
        if not bcFasta.startswith(config['metadata']):
            bcFasta = os.path.join(config['metadata'], bcFasta)
        
        if os.path.isfile(bcFasta):
            try:
                barcodes = list(SeqIO.parse(bcFasta, 'fasta'))
                if len(barcodes) == 0:
                    print(f"[WARNING] Barcode fasta file `{bcFasta}` empty or not fasta format\n", file=sys.stderr)
                if any(['_' in bc.id for bc in barcodes]):
                    print(f"[WARNING] Sequence ID(s) in barcode fasta file `{bcFasta}` contain underscore(s), which may disrupt the pipeline.\n", file=sys.stderr)
                
                # Check barcode lengths
                if 'N' in context:
                    #TODO: use the same logic as in demux and UMI_extract to determine barcode length
                    barcodeLength = context.rindex('N') - context.index('N') + 1
                    if any([len(seq.seq) != barcodeLength for seq in barcodes]):
                        print(f"[WARNING] Barcode fasta file `{bcFasta}` contains barcodes of different lengths than the context `{context}` for barcode type `{barcodeType}` in run tag `{tag}`. Demultiplexing will fail.\n", file=sys.stderr)
            except Exception as e:
                print(f"[WARNING] Failed to parse barcode fasta file `{bcFasta}`: {e}\n", file=sys.stderr)
        
        elif not barcode_data.get('generate', False):
            print(f"[WARNING] Barcode fasta file `{bcFasta}` does not exist, but is used for barcode type `{barcodeType}` in run tag `{tag}`\n", file=sys.stderr)
        
        # Validate barcode groups (both partition and label)
        for group_type in ['partition', 'label']:

            csv_key = f"{group_type}_barcode_groups_csv"

            if csv_key in config['runs'][tag]:
                try:
                    groups_df = pd.read_csv(config['runs'][tag][csv_key])
                    if 'tag' in groups_df.columns:
                        groups_df = groups_df[groups_df['tag'] == tag]
                    
                    for _, row in groups_df.iterrows():
                        group_name = row['barcode_group']
                        
                        # Check if barcode types in group are properly defined
                        bc_types = groups_df.columns.difference(['barcode_group', 'tag'])
                        for bc_type in bc_types:
                            if pd.notna(row[bc_type]):
                                if bc_type not in barcode_info:
                                    print(f"[WARNING] Barcode type `{col}` in {group_type} group `{group_name}` for run tag `{tag}` is not defined in barcode_info_csv. Demultiplexing will fail.\n", file=sys.stderr)
                                
                        # Check if barcode is used in correct group type
                        if barcodeType in row and pd.notna(row[barcodeType]):
                            # For partition groups, label_only should be False
                            # For label groups, label_only should be True
                            expected_label_only = (group_type == 'label')
                            
                            if label_only != expected_label_only:
                                if group_type == 'partition':
                                    print(f"[WARNING] `label_only` set to True for barcode type `{barcodeType}` in run tag `{tag}`, but is used for partitioning in group `{group_name}`. Demultiplexing will fail.\n", file=sys.stderr)
                                else:
                                    print(f"[WARNING] `label_only` set to False for barcode type `{barcodeType}` in run tag `{tag}`, but is used for labeling in group `{group_name}`. Demultiplexing will fail.\n", file=sys.stderr)
                            
                            # Check if barcode name exists in fasta
                            bcName = row[barcodeType]
                            if os.path.isfile(bcFasta) and (bcName not in [seq.id for seq in SeqIO.parse(bcFasta, 'fasta')]):
                                print(f"[WARNING] Barcode {bcName}, barcode type `{barcodeType}` for {group_type} group `{group_name}` for run tag `{tag}` is not present in the barcode fasta file `{bcFasta}`.\n", file=sys.stderr)
                except Exception as e:
                    print(f"[WARNING] Failed to validate {group_type} groups for tag `{tag}`: {e}\n", file=sys.stderr)
        
        # Check generate option
        if barcode_data.get('generate', False):
            num_to_generate = barcode_data['generate']
            if (num_to_generate != 'all') and not isinstance(num_to_generate, int):
                print(f"[WARNING] `generate` option for barcode type `{barcodeType}` for run tag `{tag}` is not properly defined. Must be an integer or 'all'.\n", file=sys.stderr)
            if os.path.isfile(bcFasta):
                notice = f"[NOTICE] `generate` option for barcode type `{barcodeType}` for run tag `{tag}` set to `{num_to_generate}`, but barcode fasta file `{bcFasta}` exists. Using this file for demultiplexing.\n"
                notices.append(notice)
                notice_of_notices = True
            else:
                notice = f"[NOTICE] `generate` option for barcode type `{barcodeType}` for run tag `{tag}` set to `{num_to_generate}`, and barcode fasta file `{bcFasta}` does not exist. Generating barcode fasta file containing {num_to_generate} barcodes prior to demultiplexing.\n"
                notices.append(notice)
                notice_of_notices = True
    
    # Check for duplicate contexts (same context can be used for both partition and label barcodes though)
    for bc_type, context_list in contexts.items():
        if len(set(context_list)) != len(context_list):
            print(f"[WARNING] Duplicate {bc_type} barcode contexts provided for run tag `{tag}`.\n", file=sys.stderr)


# check that tags and barcodeGroup names don't contain underscores
for tag in config['runs']:
    if '_' in tag:
        print(f"[WARNING] Run tag `{tag}` contains underscore(s), which will disrupt the pipeline. Please remove all underscores in run tag names.", file=sys.stderr)
    if 'barcodeGroups' in config['runs'][tag]:
        for bcGroup in config['runs'][tag]['barcodeGroups']:
            if '_' in bcGroup:
                print(f"[WARNING] Barcode group `{bcGroup}` for run tag `{tag}` contains underscore(s), which will disrupt the pipeline. Please remove all underscores in barcode group names.", file=sys.stderr)

# Check that 'background' partition group, if declared, is defined in all tags
if config.get('background', False):
    for tag in config['runs']:
        if 'partition_barcode_groups_csv' in config['runs'][tag]:
            try:
                groups_df = pd.read_csv(config['runs'][tag]['partition_barcode_groups_csv'])
                if 'tag' in groups_df.columns:
                    groups_df = groups_df[groups_df['tag'] == tag]
                
                group_names = groups_df['partition_bc_group'].unique()
                if config['background'] not in group_names:
                    print(f"[WARNING] `background` partition group declared in config file as {config['background']}, but this group is not defined for `{tag}`. Some pipeline rules will fail.\n", file=sys.stderr)
            except Exception:
                print(f"[WARNING] Could not verify background partition group for tag `{tag}`.\n", file=sys.stderr)
        else:
            print(f"[WARNING] `background` partition group declared in config file, but `partition_barcode_groups_csv` not supplied for run tag `{tag}`. Some pipeline rules will fail.\n", file=sys.stderr)

# add timepoints files to config dictionary in the format {'timepoints':{tag:timepointCSVfile}}.
#   This is to allow timepoint CSV files to be used once, only for the first tag that uses that file.
#   Also creates a dict to map each row to some information about that row.
#   Also checks for the following:
#       - csv file exists
#       - Referenced tags and barcodeGroups are defined in config file
#       - reference fasta files for sample/barcodeGroup combinations are the same in each row
#       - at least two timepoints are given
#       - a row only uses tags from the same sequencing run or, if it uses different sequencing runs, that a 'background'
#           barcode group is provided in the config file. This is important because background subtraction is necessary
#           for accurate rate calculations, and sequencing error can of course differ from run to run.

# Update enrichment analysis section
config['do_enrichment_analysis'] = {}

#TODO: update enrichment analysis to work for both genotypes and barcodes
if 'timepoints' in config:
    for tag in config['runs']:
        # Check for label-only barcodes
        label_only_count = 0
        if 'barcode_info_csv' in config['runs'][tag]:
            try:
                barcode_df = pd.read_csv(config['runs'][tag]['barcode_info_csv'])
                if 'tag' in barcode_df.columns:
                    barcode_df = barcode_df[barcode_df['tag'] == tag]
                label_only_count = barcode_df['label_only'].sum() if 'label_only' in barcode_df.columns else 0
            except Exception:
                pass
        
        # Enrichment analysis requires exactly 1 label-only barcode
        config['do_enrichment_analysis'][tag] = (label_only_count == 1)

    CSVpath = os.path.join(config['metadata'], config['timepoints'])
    config['timepoints'] = CSVpath

    timepoints_info_dict = {}

    # if 'timepoints' not in config: # construct timepoints dict for first tag encountered with timepoints file declared

    # if CSVpath in config['timepoints'].values():
    #     print(f"[NOTICE] Timepoints file `{CSVpath}` is used for multiple tags. Files will be named only using the first tag that uses this file.\n", file=sys.stderr)
    # else:
    #     config['timepoints'][tag] = CSVpath

    # decide whether to do barcode enrichment analysis

    if os.path.isfile(CSVpath):
        timepointsCSV = pd.read_csv(CSVpath, index_col=0, header=1)
        topRow = [x for x in pd.read_csv(CSVpath).columns if 'Unnamed: ' not in x]

        if len(topRow) > 1:
            print(f"[NOTICE] More than one cell is filled in the top row of timepoint CSV file {CSVpath}. Only the first cell in this row will be used for labeling timepoint-based plots.\n", file=sys.stderr)
            units = topRow[0]
        elif len(topRow) == 0: 
            print(f"[NOTICE] No time unit provided in top row of timepoint CSV file {CSVpath}. Default 'generations' will be used.\n", file=sys.stderr)
            units = 'generations'
        else:
            units = topRow[0]

        if len(timepointsCSV.columns) <= 1:
            print(f"[WARNING] Timepoints .CSV file `{CSVpath}` does not have at least two timepoints. Timepoint-based snakemake rules will fail.\n", file=sys.stderr)
        else:
            rowIndex = 3    # start at 3 because first two rows are ignored with pd.read_csv call, and errors/warnings will use 1-indexing
            uniqueSamples = list(timepointsCSV.index.unique())
            replicates = False
            if len(uniqueSamples) != len(timepointsCSV.index): # assign replicate designations for samples with replicates (more than one row with the same sample name)
                replicates = True

            for sampleGroup in uniqueSamples:

                replicateIndex = 0

                if sampleGroup in config['runs']:
                    errors.append(f"[ERROR] Timepoint sample name used for row {rowIndex} (1-indexed) of timepoints .CSV file `{CSVpath}`, `{sampleGroup}` is also used as a tag. This is forbidden.\n")
                if '_' in sampleGroup:
                    errors.append(f"[ERROR] Underscore used in timepoint sample name `{sampleGroup}` in row {rowIndex} (1-indexed) of timepoints .CSV file `{CSVpath}`. This is forbidden.\n")
                # if (sampleGroup in config['timepointsInfo']) or (sampleGroup in [key.split('-rep')[0] for key in config['timepointsInfo']]):
                #     print(f'[NOTICE] Sample name `{sampleGroup}` is used in timepoints file `{CSVpath}` as well as in another timepoints file. This is allowed, but note that plots in `plots/timepoints` will only be based on the first timepoints file that uses this sample name.\n', file=sys.stderr)
                #     continue

                for _, row in timepointsCSV.loc[timepointsCSV.index==sampleGroup].iterrows():
                    replicateIndex += 1
                    if replicates:
                        sample = sampleGroup + f"-rep{str(replicateIndex)}"
                    else:
                        sample = sampleGroup

                    firstTP = timepointsCSV.columns[0]
                    i = 0
                    while pd.isnull(row[firstTP]): # not all samples need to have the same first timepoint
                        i+=1
                        firstTP = timepointsCSV.columns[i]
                    firstTag = str(row[firstTP]).split('_')[0]
                    if firstTag in config['runs']:
                        firstTagRefFasta = config['runs'][firstTag]['reference']
                        if config['do_NT_mutation_analysis'].get(firstTag, False): # only necessary if NT analysis is being run
                            timepoints_info_dict[sample] = {'tag':firstTag, 'reference':firstTagRefFasta, 'units':units, 'tag_barcode_tp':{}}
                        firstTagRefSeq = str(list(SeqIO.parse(firstTagRefFasta, 'fasta'))[0].seq).upper()
                        firstTagRunname = config['runs'][firstTag].get('runname', '')
                    else:
                        errors.append(f"[ERROR] Tag referenced in row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{firstTag}` is not defined in config file. Check timepoints csv file and config file for errors.\n")
                    for tp in timepointsCSV.columns:
                        if str(row[tp]) == 'nan':
                            continue
                        tag_barcodeGroup = row[tp].split('_')
                        if len(tag_barcodeGroup) != 2:
                            errors.append(f"[ERROR] Timepoint sample `{sample}` in row {rowIndex} of timepoints .CSV file `{CSVpath}` does not have a valid tag_barcodeGroup format. Must be in the format `tag_barcodeGroup`.\n")
                            continue
                        tag, barcodeGroup = tag_barcodeGroup
                        if tag in config['runs']:
                            if firstTag in config['runs']:
                                tagRefSeq = str(list(SeqIO.parse(config['runs'][tag]['reference'], 'fasta'))[0].seq).upper()
                                if tagRefSeq != firstTagRefSeq:
                                    print(f"[WARNING] In row {rowIndex} of timepoints .CSV file `{CSVpath}`, samples `{row[firstTP]}` and `{row[tp]}` use different reference sequences. Analysis may be unreliable.\n", file=sys.stderr)
                                tagRunname = config['runs'][tag].get('runname', '')
                                if tagRunname != firstTagRunname and 'background' not in config:
                                    print(f"[WARNING] In row {rowIndex} of timepoints .CSV file `{CSVpath}`, samples `{row[firstTP]}` and `{row[tp]}` use different runnames, but a background barcodeGroup is not provided. Analysis may be unreliable.\n", file=sys.stderr)
                            if config['do_NT_mutation_analysis'].get(firstTag, False): # only necessary if NT analysis is being run
                                timepoints_info_dict[sample]['tag_barcode_tp'][(tag, barcodeGroup)] = tp
                        else:
                            errors.append(f"[ERROR] Tag referenced in row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{tag}` is not defined in config file. Check timepoints csv file and config file for errors\n")
            config['timepointsInfo'] = timepoints_info_dict
    else:
        print(f"[WARNING] Timepoints .CSV file for run tag `{tag}`, `{CSVpath}` does not exist.\n", file=sys.stderr)
if len(errors) > 0:
    for err in errors:
        print(err, file=sys.stderr)
    raise RuntimeError("Critical errors found. See above.\n")

if (rename_csv := config.get('rename', False)):
    rename_path = os.path.join(config['metadata'], rename_csv)
    if os.path.isfile(rename_path):
        rename_df = pd.read_csv(rename_path)
        config['rename'] = dict(zip(rename_df.iloc[:, 1], rename_df.iloc[:, 0]))
    else:
        print(f"[WARNING] Rename .CSV file `{rename_csv}` does not exist.\n", file=sys.stderr)

if notice_of_notices:
    print("[NOTICE] Some notices not printed to terminal. Please see log file for more information.\n", file=sys.stderr)

with open(os.path.join(config['metadata'], 'config_final.yaml'), 'w') as outfile:
    yaml.dump(config, outfile, default_flow_style=False)

# check for required options
missing = []
for option in required:
    if option not in config:
        missing.append(option)
if len(missing) > 0:
    text = [f"`{o}`" for o in missing]
    print(f"[WARNING] Required option(s) missing from the config file: {', '.join(text)}. Please add these options to the config file. See example_working_directory/config.yaml for example.\n", file=sys.stderr)

# # include modules
include : "rules/get_seqs.smk"
include : "rules/consensus.smk"
include : "rules/analysis.smk"
include : "rules/target_files.smk"
include : "rules/clean.smk"

# error and success handler
def print_log(status='SUCCESS'):
    os.makedirs('log', exist_ok=True)
    now = datetime.now()
    log_name = os.path.join('log', now.strftime('%Y%m%d_%H_%M_%S_%f.maple.log'))
    end_files = get_dir_files(workflow.workdir_init)
    sm_log = None
    for h in logging.getLogger().handlers:
        if isinstance(h, logging.FileHandler):
            sm_log = h.baseFilename
            break
    if sm_log:
        sm_log = os.path.relpath(sm_log)
    else:
        sm_log = '<unknown>'
    with open(log_name, 'w') as fp:
        print('Log file for maple version {tag}'.format(tag=maple_tag), file=fp)
        print("Workflow begin: {}".format(start_time.strftime('%d.%m.%Y %H:%M:%S')), file=fp)
        print("Workflow end:   {}".format(now.strftime('%d.%m.%Y %H:%M:%S')), file=fp)
        print('Command: {}'.format(' '.join(sys.argv)), file=fp)
        print('', file=fp)
        print("Status: {}".format(status), file=fp)
        print('', file=fp)
        print("Working directory: {}".format(workflow.workdir_init), file=fp)
        print("Log file: {}".format(log_name), file=fp)
        print("Snakemake log file: {}".format(os.path.relpath(sm_log)), file=fp)
        print('', file=fp)
        print("maple config:", file=fp)
        print('-----------------------------------', file=fp)
        print(yaml.dump({key:value for key, value in config.items()
            if not (isinstance(value, dict) or isinstance(value, list))}, indent=2, sort_keys=True), file=fp)
        print("Environment config", file=fp)
        print('-----------------------------------', file=fp)
        print(yaml.dump({key:value for key, value in config.items()
            if isinstance(value, dict) or isinstance(value, list)}, indent=2, sort_keys=True), file=fp)
        print("File system changes:", file=fp)
        print('-----------------------------------', file=fp)
        print("New files:", file=fp)
        print('\n'.join(sorted([f for f in end_files.difference(start_files)])), file=fp)
        print("Deleted files:", file=fp)
        print('\n'.join(sorted([f for f in start_files.difference(end_files)])), file=fp)
        if len(notices) > 0:
            print('Notices:', file=fp)
            print('-----------------------------------', file=fp)
            print('\n'.join(notices), file=fp)
    return log_name

onsuccess:
    log_name = print_log(status='SUCCESS')
    print("""
maple completed successfully.
the log file was written to {}.""".format(log_name), file=sys.stderr)


onerror:
    log_name = print_log(status='ERROR')
    print("""
maple exited with an error.
the log file was written to {}.

please visit the github at
    https://github.com/gordonrix/maple
to make sure everything is configured correctly.

""".format(log_name), file=sys.stderr)
