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
import pandas as pd

start_time = datetime.now()


# snakemake config
min_version("5.5.2")
configfile: "config.yaml"

# filter output depending on run mode
def print_(*args, **kwargs):
    if workflow.mode == snakemake.common.Mode.default:
        print(*args, **kwargs)

print('')


# get pipeline version # update for maple
def get_tag():
    try:
        cmd = 'git describe --tags'
        version = subprocess.check_output(cmd.split(), cwd=os.path.dirname(workflow.snakefile)).decode().strip()
    except subprocess.CalledProcessError:
        print_('[WARNING] Unable to get version from git tags.', file=sys.stderr)
        version = '-'
    try:
        cmd = 'git rev-parse --abbrev-ref HEAD'
        branch = subprocess.check_output(cmd.split(), cwd=os.path.dirname(workflow.snakefile)).decode().strip()
    except subprocess.CalledProcessError:
        print_('[WARNING] Unable to get branch from git. Pulling development.', file=sys.stderr)
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


# append username to shadow prefix if not present
if hasattr(workflow, "shadow_prefix") and workflow.shadow_prefix:
    shadow_prefix = workflow.shadow_prefix
    if not os.environ['USER'] in shadow_prefix:
        shadow_prefix = os.path.join(shadow_prefix, os.environ['USER'])
        print_("[INFO] Shadow prefix is changed from {p1} to {p2} to be user-specific".format(
            p1=workflow.shadow_prefix, p2=shadow_prefix), file=sys.stderr)
    workflow.shadow_prefix = shadow_prefix


# parse pipeline environment
maple_env = {}
with open(os.path.join(os.path.dirname(workflow.snakefile), "env.yaml"), 'r') as fp:
    maple_env = yaml.safe_load(fp)


# verify given binaries
if 'bin' in maple_env:
    if not 'bin' in config:
        config['bin'] = {}
    if not 'bin_singularity' in config:
        config['bin_singularity'] = {}
    for name, loc in maple_env['bin'].items():
        loc_sys = None
        loc_singularity = os.path.basename(loc)
        if os.path.isfile(loc):
            # absolute path is given
            loc_sys = loc
        elif os.path.isfile(os.path.join(os.path.dirname(sys.executable), loc)):
            # executable in python installation/virtual environment
            loc_sys = os.path.join(os.path.dirname(sys.executable), loc)
        else:
            # scan the PATH if we find the executable
            for path in os.environ["PATH"].split(os.pathsep):
                f = os.path.join(path, os.path.basename(loc))
                if os.path.isfile(f):
                    loc_sys = f
                    break
        # save executable path depending on singularity usage
        if hasattr(workflow, 'use_singularity') and workflow.use_singularity:
            # within singularity everything is accessible through /bin and /usr/bin
            config['bin_singularity'][name] = loc_singularity
            if loc_sys:
                config['bin'][name] = loc_sys
            #else:
            #    print_("[WARNING] {name} not found as {loc} and is only available in singularity rules.".format(
            #        name=name, loc=loc), file=sys.stderr)
        else:
            # singularity rules use system wide executables
            if loc_sys:
                config['bin_singularity'][name] = loc_sys
                config['bin'][name] = loc_sys
            else:
                print_("[WARNING] {name} not found as {loc} and is not available in the workflow.".format(
                    name=name, loc=loc), file=sys.stderr)
else:
    raise RuntimeError("[ERROR] No binaries in environment configuration\n", file=sys.stderr)


# Runtime scaling of data depending tools
if 'runtime' in maple_env:
    config['runtime'] = {}
    for key, value in maple_env['runtime'].items():
        config['runtime'][key] = value
else:
    raise RuntimeError("[ERROR] No runtime scalings in environment configuration\n", file=sys.stderr)


# memory scaling of data depending tools
if 'memory' in maple_env:
    config['memory'] = {}
    for key, value in maple_env['memory'].items():
        config['memory'][key] = tuple(value)
else:
    raise RuntimeError("[ERROR] No memory scalings in environment configuration.\n", file=sys.stderr)



# locations of helper scripts in rules/utils
if not 'sbin' in config:
    config['sbin'] = {}
if not 'sbin_singularity' in config:
    config['sbin_singularity'] = {}
for s in [s for s in os.listdir(os.path.join(os.path.dirname(workflow.snakefile), 'rules/utils/')) if
        os.path.isfile(os.path.join(os.path.dirname(workflow.snakefile), 'rules/utils', s))]:
    if s.startswith('__') or s.startswith('.'):
        continue
    config['sbin'][s] = os.path.join(os.path.dirname(workflow.snakefile), 'rules/utils', s)
    if hasattr(workflow, 'use_singularity') and workflow.use_singularity:
        config['sbin_singularity'][s] = os.path.join('/app/rules/utils', s)
    else:
        config['sbin_singularity'][s] = config['sbin'][s]


# helper of submodules are called relative to the pipeline base directory
config['sbin']['base'] = os.path.join(os.path.dirname(workflow.snakefile))
if hasattr(workflow, 'use_singularity') and workflow.use_singularity:
    config['sbin_singularity']['base'] = '/app'
else:
    config['sbin_singularity']['base'] = config['sbin']['base']


# find the python executable
# Python executable of the workflow
config['bin']['python'] = sys.executable

# In the container we just use python3
if hasattr(workflow, 'use_singularity') and workflow.use_singularity:
	config['bin_singularity']['python'] = 'python3'
else:
	config['bin_singularity']['python'] = sys.executable

errors = []

# list of required options. Will be added to if certain tools are used
required = ['sequences_dir', 'fastq_dir', 'references_directory', 'threads_alignment', 'threads_samtools', 'threads_demux', 'NGmerge_flags', 'nanopore', 'nanoplot', 'nanoplot_flags', 'alignment_samtools_flags', 'alignment_minimap2_flags', 'mutation_analysis_quality_score_minimum', 'sequence_length_threshold', 'highest_abundance_genotypes', 'mutations_frequencies_raw', 'unique_genotypes_count_threshold', 'runs']
if ('minknowDir' in config) and ('sequences_dir' not in config):
    print('[WARNING] "sequences_dir" is a new required config variable that replaces "minknowDir", please replaced "minknowDir" with "sequences_dir" in your config file"')

config['merge_paired_end'] = {}
for tag in config['runs']:
    config['merge_paired_end'][tag] = False
    if any([x in config['runs'][tag] for x in ['fwdReads', 'rvsReads']]):
        config['merge_paired_end'][tag] = True

runs_to_import = []
# check for sequences
for tag in config['runs']:
    sequences = os.path.join('sequences', tag+'.fastq.gz')

    if not os.path.exists(sequences):

        if config['merge_paired_end'][tag] == True:
            if 'fwdReads' not in config['runs'][tag] or 'rvsReads' not in config['runs'][tag]:
                print_(f"[WARNING] Both forward and reverse reads files not provided for {tag} with keyword `fwdReads` and `rvsReads`.\n", file=sys.stderr)
            fwd = os.path.join('sequences', 'paired', config['runs'][tag]['fwdReads'])
            rvs = os.path.join('sequences', 'paired', config['runs'][tag]['rvsReads'])
            if not all((os.path.exists(fwd), os.path.exists(rvs))):
                print_(f"[WARNING] One or both forward/reverse reads files provided for {tag}, ({fwd}, {rvs}) do not exist.\n", file=sys.stderr)

        elif 'runname' not in config['runs'][tag]:
            print_(f"[WARNING] Neither paired end reads nor runname director(y/ies) provided for tag `{tag}`, but sequences file `{sequences}` not found.\n", file=sys.stderr)

        else:
            for runname in config['runs'][tag]['runname']:
                batch = os.path.join('sequences', 'batches', runname)
                if not os.path.exists(batch):
                    runs_to_import.append(runname)
# Check minknow directory
if runs_to_import != []:
    for runname in runs_to_import:
        if not os.path.isdir(config['sequences_dir']):
            print_(f"[WARNING] May need to import runname `{runname}`, but the provided minknow directory, `{config['sequences_dir']}`, does not exist.\n", file=sys.stderr)
        elif all([runname not in dirs for _, dirs, _ in os.walk(config['sequences_dir'].rstrip('/'))]):
            print_(f"[WARNING] May need to import runname `{runname}`, but this could not be located in any directory tree under `{config['sequences_dir']}`.\n", file=sys.stderr)


# check reference sequences
refSeqFastaFiles = []   # list files for all tags then check so files not being checked multiple times
config['do_NT_mutation_analysis'] = {}     # dictionaries to determine if NT/AA analysis should be performed based on how many sequences are present in the ref fasta file
config['do_AA_mutation_analysis'] = {}
for tag in config['runs']:
    if 'reference' not in config['runs'][tag]:
        errors.append(f"[ERROR] No reference file provided for tag `{tag}")
    refName = config['runs'][tag]['reference']
    refFullPath = os.path.join(config['references_directory'], config['runs'][tag]['reference'])
    if not (refName.endswith('fasta') or refName.endswith('.fa')):
        print_(f'[WARNING] Reference .fasta file for {tag} does not end with `.fasta` or `.fa` (given path: {refFullPath}).', file=sys.stderr)
    config['runs'][tag]['reference'] = refFullPath
    if not os.path.isfile(refFullPath):
        print_(f'[ERROR] Reference .fasta file for {tag} (given path: {refFullPath}) not found.', file=sys.stderr)
    refFastaPrefix = refName.split('.f')[0]
    alnRefFullPath = os.path.join(config['references_directory'], '.' + refFastaPrefix + '_aln.fasta')
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
        print_(f'[WARNING] AA_muts_of_interest provided for run tag `{tag}`, but no protein seq provided for this tag. AA muts of interest will not be evaluated for this tag.', file=sys.stderr)

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
            if nt not in list("ATGC"):
                errors.append(f"[ERROR] Character {nt} at position {i} in reference sequence `{proteinSeq.id}` of reference file `{refFasta}` is not a canonical nucleotide\n")

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
        print_(f'Alignment reference .fasta file not found or is different from original reference .fasta file. Generating {alnFasta} from {refFasta}.\n', file=sys.stderr)
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
        splintFasta = os.path.join(config['references_directory'], f".{tag}_splint.fasta")
        makeSplintFasta = True
        if os.path.isfile(splintFasta):
            first_record = next(SeqIO.parse(splintFasta, 'fasta'))
            if (first_record.id == 'splint') and (str(first_record.seq).upper() == splintSeq.upper()):
                makeSplintFasta = False
            else:
                print_(f"[NOTICE] Splint sequence provided in config file for run tag `{tag}` has changed. Updating splint fasta file `{splintFasta}`.\n", file=sys.stderr)
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
            print_(f"[NOTICE] Raw data / alignment sequence / UMI context / splint combination used more than once. Using the consensus .fasta file of tag `{consensusRecipeDict[consensusRecipe]}` for tag `{tag}` to reduce computation time and storage requirements.\n", file=sys.stderr)
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
    if 'barcodeInfo' not in config['runs'][tag]:
        print_(f"[NOTICE] No keyword 'barcodeInfo' provided for {tag}, will not perform demultiplexing for this tag.\n", file=sys.stderr)
        config['do_demux'][tag] = False
        continue
    else:
        config['do_demux'][tag] = True
    if len(config['runs'][tag]['barcodeInfo']) == 0:
        print_(f"[WARNING] `barcodeInfo` for run tag `{tag}` does not contain any barcode types. Demultiplexing will fail.\n", file=sys.stderr)
    if 'barcodeGroups' in config['runs'][tag]:
        # add barcodeGroups to tag as a dict if declared as a csv file
        if type(config['runs'][tag]['barcodeGroups']) == str:
            CSVpath = os.path.join(config['references_directory'], config['runs'][tag]['barcodeGroups'])
            if os.path.isfile(CSVpath):
                barcodeGroupsCSV = pd.read_csv(CSVpath, index_col=False, header=1, dtype=str)
                barcodeGroupsCSV = barcodeGroupsCSV.set_index(barcodeGroupsCSV.columns[0]) #can't do this in one line because of a pd.read_csv bug that doesn't allow index to be string
                if barcodeGroupsCSV.index.name == 'tag':
                    barcodeGroupsCSV = barcodeGroupsCSV.loc[barcodeGroupsCSV.index==tag].set_index('barcodeGroup')
                if any([c[:9]=='Unnamed: ' for c in barcodeGroupsCSV.columns]):
                    print_(f"[WARNING] Barcode type beginning with 'Unnamed: ' detected for tag {tag} in barcodeGroups csv {CSVpath}. This usually results from erroneous whitespace characters. Demultiplexing may fail.\n", file=sys.stderr)
                config['runs'][tag]['barcodeGroups'] = barcodeGroupsCSV.to_dict('index')
            else:
                print_(f"[NOTICE] String provided for `barcodeGroups` in run tag `{tag}`, but file path `{CSVpath}` does not exist. Will use barcode combinations to name demultiplexed files.", file=sys.stderr)
        if len(config['runs'][tag]['barcodeGroups']) == 0:
            print_(f"[NOTICE] No barcode groups provided for run tag `{tag}`. Outputs will be named as concatemerized barcode names.\n", file=sys.stderr)
    else:
        print_(f"[NOTICE] `barcodeInfo` supplied but `barcodeGroups` not supplied as dict or .CSV file for run tag `{tag}`. Will use barcode combinations to name demultiplexed files.\n", file=sys.stderr)
    refFasta = config['runs'][tag]['reference']
    alignmentSeq = list(SeqIO.parse(refFasta, 'fasta'))[0]
    contexts = []
    for barcodeType in config['runs'][tag]['barcodeInfo']:
        for requiredKey in ['context', 'fasta', 'reverseComplement']:
            if requiredKey not in config['runs'][tag]['barcodeInfo'][barcodeType]:
                print_(f"[WARNING] Tag `{tag}` barcode type `{barcodeType}` does not contain the required key `{requiredKey}`.\n", file=sys.stderr)
        c = config['runs'][tag]['barcodeInfo'][barcodeType].get('context', False)
        if c: contexts.append(c)
        config['runs'][tag]['barcodeInfo'][barcodeType]['context'] = config['runs'][tag]['barcodeInfo'][barcodeType]['context'].upper()
        occurences = str(alignmentSeq.seq).upper().count(config['runs'][tag]['barcodeInfo'][barcodeType]['context'].upper())
        if occurences == 0:
            print_(f"[WARNING] Barcode type `{barcodeType}` context `{config['runs'][tag]['barcodeInfo'][barcodeType]['context']}` not found in reference `{alignmentSeq.id}` in fasta `{refFasta}`\n", file=sys.stderr)
        elif occurences > 1:
            print_(f"[WARNING] Barcode type `{barcodeType}` context `{config['runs'][tag]['barcodeInfo'][barcodeType]['context']}` found more than once in reference `{alignmentSeq.id}` in fasta `{refFasta}`\n", file=sys.stderr)
        bcFasta = os.path.join(config['references_directory'], config['runs'][tag]['barcodeInfo'][barcodeType]['fasta'])
        config['runs'][tag]['barcodeInfo'][barcodeType]['fasta'] = bcFasta
        if os.path.isfile(bcFasta):
            if len(list(SeqIO.parse(bcFasta, 'fasta'))) == 0:
                print_(f"[WARNING] Barcode fasta file `{bcFasta}` empty or not fasta format\n\n", file=sys.stderr)
            if any(['_' in bc.id for bc in list(SeqIO.parse(bcFasta, 'fasta'))]):
                print_(f"[WARNING] Sequence ID(s) in barcode fasta file `{bcFasta}` contain underscore(s), which may disrupt the pipeline. Please remove all underscores in sequence IDs.", file=sys.stderr)
            if type(config['runs'][tag]['barcodeInfo'][barcodeType]['reverseComplement'])!=bool:
                print_(f"[WARNING] Tag `{tag}`, barcode type `{barcodeType}` reverseComplement keyword must be set as True or False\n\n", file=sys.stderr)
        elif config['runs'][tag]['barcodeInfo'][barcodeType].get('generate', False) == False:
            print_(f"[WARNING] Barcode fasta file `{bcFasta}` does not exist, but is used for barcode type `{barcodeType}` in run tag `{tag}`\n", file=sys.stderr)
        if 'barcodeGroups' in config['runs'][tag]:
            for bcGroup in config['runs'][tag]['barcodeGroups']:
                for bcType in config['runs'][tag]['barcodeGroups'][bcGroup]:
                    if bcType not in config['runs'][tag]['barcodeInfo']:
                        print_(f"[WARNING] Barcode type `{bcType}` in barcode group `{bcGroup}` for run tag `{tag}` is not defined in 'barcodeInfo'. Demultiplexing will fail.\n", file=sys.stderr)
                if config['runs'][tag]['barcodeInfo'][barcodeType].get('noSplit', False) == True:
                    for bcType in config['runs'][tag]['barcodeGroups'][bcGroup]:
                        if bcType == barcodeType:
                            print_(f"[WARNING] `noSplit` set to True for barcode type `{barcodeType}` in run tag `{tag}`, but is used for naming in barcode group `{bcGroup}`. Demultiplexing will fail.\n", file=sys.stderr)
                elif config['runs'][tag]['barcodeInfo'][barcodeType].get('noSplit', False) == False:
                    if os.path.isfile(bcFasta) and (config['runs'][tag]['barcodeGroups'][bcGroup][barcodeType] not in [seq.id for seq in list(SeqIO.parse(bcFasta, 'fasta'))]):
                        print_(f"[WARNING] Barcode type `{barcodeType}` in barcode group `{bcGroup}` for run tag `{tag}` is not present in the barcode fasta file `{config['runs'][tag]['barcodeInfo'][barcodeType]['fasta']}` set for this tag.\n", file=sys.stderr)
        if 'generate' in config['runs'][tag]['barcodeInfo'][barcodeType]:
            numToGenerate = config['runs'][tag]['barcodeInfo'][barcodeType]['generate']
            if (numToGenerate != 'all') and type(numToGenerate) != int:
                print_(f"[WARNING] `generate` option for barcode type `{barcodeType}` for run tag `{tag}` is not properly defined. Must be an integer or 'all'.\n", file=sys.stderr)
            if os.path.isfile(bcFasta):
                print_(f"[NOTICE] `generate` option for barcode type `{barcodeType}` for run tag `{tag}` set to `{numToGenerate}`, but barcode fasta file `{config['runs'][tag]['barcodeInfo'][barcodeType]['fasta']}` exists. Using this file for demultiplexing.\n", file=sys.stderr)
            else:
                print_(f"[NOTICE] `generate` option for barcode type `{barcodeType}` for run tag `{tag}` set to `{numToGenerate}`, and barcode fasta file `{config['runs'][tag]['barcodeInfo'][barcodeType]['fasta']}` does not exist. Generating barcode fasta file containing {numToGenerate} barcodes prior to demultiplexing.\n", file=sys.stderr)
        if len(set(contexts)) != len(contexts):
            print_(f"[WARNING] Duplicate barcode contexts provided for run tag `{tag}`.\n", file=sys.stderr)

# check that tags and barcodeGroup names don't contain underscores
for tag in config['runs']:
    if '_' in tag:
        print_(f"[WARNING] Run tag `{tag}` contains underscore(s), which will disrupt the pipeline. Please remove all underscores in run tag names.", file=sys.stderr)
    if 'barcodeGroups' in config['runs'][tag]:
        for bcGroup in config['runs'][tag]['barcodeGroups']:
            if '_' in bcGroup:
                print_(f"[WARNING] Barcode group `{bcGroup}` for run tag `{tag}` contains underscore(s), which will disrupt the pipeline. Please remove all underscores in barcode group names.", file=sys.stderr)

# check that 'background' barcodeGroup, if declared, is defined in all tags:
if config.get('background', False):
    for tag in config['runs']:
        if 'barcodeGroups' in config['runs'][tag]:
            if config['background'] not in config['runs'][tag]['barcodeGroups']:
                print_(f"[WARNING] `background` barcodeGroup declared in config file as {config['background']}, but this barcodeGroup is not defined for `{tag}`. Some pipeline rules will fail.\n", file=sys.stderr)
        else:
            print_(f"[WARNING] `background` barcodeGroup declared in config file, but `barcodeGroups` not supplied as dict or .CSV file for run tag `{tag}`. Some pipeline rules will fail.\n", file=sys.stderr)

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
config['do_enrichment_analysis'] = {}
for tag in config['runs']:
    if 'timepoints' in config['runs'][tag]:
        CSVpath = os.path.join(config['references_directory'], config['runs'][tag]['timepoints'])
        if 'timepoints' not in config: # construct timepoints dict for first tag encountered with timepoints file declared
            config['timepoints'] = {}
            config['timepointsInfo'] = {}
        if CSVpath not in config['timepoints'].values():
            config['timepoints'][tag] = CSVpath

        # decide whether to do barcode enrichment analysis
        no_split_bcs = [bc for bc in config['runs'][tag]['barcodeInfo'] if config['runs'][tag]['barcodeInfo'][bc].get('noSplit', False) == True]
        # if a tag is defined with exactly 1 nosplit barcode, and a timepoints file, then enrichment analysis will be performed on that tag
        if len(no_split_bcs) == 1:
            config['do_enrichment_analysis'][tag] = True

        if os.path.isfile(CSVpath):
            timepointsCSV = pd.read_csv(CSVpath, index_col=0, header=1)
            topRow = [x for x in pd.read_csv(CSVpath).columns if 'Unnamed: ' not in x]

            if len(topRow) > 1:
                print_(f"[NOTICE] More than one cell is filled in the top row of timepoint CSV file {CSVpath}. Only the first cell in this row will be used for labeling timepoint-based plots.\n", file=sys.stderr)
            elif len(topRow) == 0: 
                print_(f"[NOTICE] No time unit provided in top row of timepoint CSV file {CSVpath}. Default 'generations' will be used.\n", file=sys.stderr)
            else:
                units = topRow[0]

            if len(timepointsCSV.columns) <= 1:
                print_(f"[WARNING] Timepoints .CSV file for run tag `{tag}`, `{CSVpath}` does not have at least two timepoints. Timepoint-based snakemake rules will fail.\n", file=sys.stderr)
            else:
                rowIndex = 3    # start at 3 because first two rows are ignored with pd.read_csv call, and errors/warnings will use 1-indexing
                uniqueSamples = list(timepointsCSV.index.unique())
                replicates = False
                if len(uniqueSamples) != len(timepointsCSV.index): # assign replicate designations for samples with replicates (more than one row with the same sample name)
                    replicates = True

                for sampleGroup in uniqueSamples:

                    replicateIndex = 0

                    if sampleGroup in config['runs']:
                        errors.append(f"[ERROR] Timepoint sample name used for row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{sampleGroup}` is also used as a tag. This is forbidden.\n")
                    if '_' in sampleGroup:
                        errors.append(f"[ERROR] Underscore used in timepoint sample name `{sampleGroup}` in row {rowIndex} of timepoints .CSV file `{CSVpath}`. This is forbidden.\n")

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
                                config['timepointsInfo'][sample] = {'tag':firstTag, 'reference':firstTagRefFasta, 'units':units, 'tag_barcode_tp':{}}
                            firstTagRefSeq = str(list(SeqIO.parse(firstTagRefFasta, 'fasta'))[0].seq).upper()
                            firstTagRunname = config['runs'][firstTag].get('runname', '')
                        else:
                            errors.append(f"[ERROR] Tag referenced in row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{firstTag}` is not defined in config file. Check timepoints csv file and config file for errors.\n")
                        for tp in timepointsCSV.columns:
                            if str(row[tp]) == 'nan':
                                continue
                            tag, barcodeGroup = row[tp].split('_')
                            if tag in config['runs']:
                                if firstTag in config['runs']:
                                    tagRefSeq = str(list(SeqIO.parse(config['runs'][tag]['reference'], 'fasta'))[0].seq).upper()
                                    if tagRefSeq != firstTagRefSeq:
                                        print_(f"[WARNING] In row {rowIndex} of timepoints .CSV file `{CSVpath}`, samples `{row[firstTP]}` and `{row[tp]}` use different reference sequences. Analysis may be unreliable.\n", file=sys.stderr)
                                    tagRunname = config['runs'][tag].get('runname', '')
                                    if tagRunname != firstTagRunname and 'background' not in config:
                                        print_(f"[WARNING] In row {rowIndex} of timepoints .CSV file `{CSVpath}`, samples `{row[firstTP]}` and `{row[tp]}` use different runnames, but a background barcodeGroup is not provided. Analysis may be unreliable.\n", file=sys.stderr)
                                if config['do_NT_mutation_analysis'].get(firstTag, False): # only necessary if NT analysis is being run
                                    config['timepointsInfo'][sample]['tag_barcode_tp'][(tag, barcodeGroup)] = tp
                            else:
                                errors.append(f"[ERROR] Tag referenced in row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{tag}` is not defined in config file. Check timepoints csv file and config file for errors\n")
                     
        else:
            print_(f"[WARNING] Timepoints .CSV file for run tag `{tag}`, `{CSVpath}` does not exist.\n", file=sys.stderr)
if len(errors) > 0:
    for err in errors:
        print_(err, file=sys.stderr)
    raise RuntimeError("Critical errors found. See above.\n")

# check for required options
missing = []
for option in required:
    if option not in config:
        missing.append(option)
if len(missing) > 0:
    text = [f"`{o}`" for o in missing]
    print_(f"[WARNING] Required option(s) missing from the config file: {', '.join(text)}. Please add these options to the config file. See example_working_directory/config.yaml for example.\n", file=sys.stderr)

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
        print("Snakemake log file: {}".format(os.path.relpath(logger.logfile)), file=fp)
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
    return log_name

onsuccess:
    if workflow.mode == snakemake.common.Mode.default:
        log_name = print_log(status='SUCCESS')
        print("""
maple completed successfully.
the log file was written to {}.""".format(log_name), file=sys.stderr)


onerror:
    if workflow.mode == snakemake.common.Mode.default:
        log_name = print_log(status='ERROR')
        print("""
maple exited with an error.
the log file was written to {}.

please visit the github at
    https://github.com/gordonrix/maple
to make sure everything is configured correctly.

""".format(log_name), file=sys.stderr)
