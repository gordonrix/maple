# \HEADER\-------------------------------------------------------------------------
#
#  CONTENTS      : Snakemake nanopore data pipeline
#
#  DESCRIPTION   : none
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#

# imports
import os, sys, collections
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
    raise RuntimeError("[ERROR] No binaries in environment configuration.")


# Runtime scaling of data depending tools
if 'runtime' in maple_env:
    config['runtime'] = {}
    for key, value in maple_env['runtime'].items():
        config['runtime'][key] = value
else:
    raise RuntimeError("[ERROR] No runtime scalings in environment configuration.")


# memory scaling of data depending tools
if 'memory' in maple_env:
    config['memory'] = {}
    for key, value in maple_env['memory'].items():
        config['memory'][key] = tuple(value)
else:
    raise RuntimeError("[ERROR] No memory scalings in environment configuration.")



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


if config['do_basecalling'] and config['merge_paired_end']:
    raise RuntimeError("[ERROR] `do_basecalling` and `merge_paired_end` cannot both be True. Set one of these to False.")

# check for required options
required = ['storage_data_raw', 'fast5_dir', 'storage_runname', 'do_basecalling', 'basecalling_guppy_config', 'basecalling_guppy_qscore_filter', 'basecalling_guppy_flags', 'porechop', 'medaka_model', 'medaka_conensus_flags', 'medaka_stitch_flags', 'references_directory', 'threads_basecalling', 'threads_porechop', 'threads_medaka', 'threads_alignment', 'threads_samtools', 'threads_demux', 'merge_paired_end', 'NGmerge_flags', 'nanopore', 'nanoplot_flags', 'UMI_consensus', 'UMI_mismatches', 'UMI_consensus_minimum', 'UMI_consensus_maximum', 'alignment_samtools_flags', 'alignment_minimap2_flags', 'demux', 'demux_screen_failures', 'demux_threshold', 'mutation_analysis_quality_score_minimum', 'sequence_length_threshold', 'do_AA_analysis', 'auto_detect_longest_ORF', 'highest_abundance_genotypes', 'mutations_frequencies_raw', 'analyze_seqs_w_frameshift_indels', 'unique_genotypes_count_threshold', 'NT_distribution_plot_x_max', 'AA_distribution_plot_x_max', 'runs']
missing = []
for option in required:
    if option not in config:
        missing.append(option)
if len(missing) > 0:
    text = [f"`{o}`" for o in missing]
    print_(f"[WARNING] Required option(s) missing from the config file: {', '.join(text)}. Please add these options to the config file. See example_working_directory/config.yaml for example.")

# check raw data archive
if config['do_basecalling']:
    if not os.path.exists(config['storage_data_raw']):
        raise RuntimeError("[ERROR] Raw data archive not found.")
    config['storage_data_raw'] = config['storage_data_raw'].rstrip('/')
    for tag in config['runs']:
        for runname in config['runs'][tag]['runname']:
            loc = os.path.join(config['storage_data_raw'], runname)
            if not os.path.exists(loc):
                print_("[WARNING] {runname} not found at {loc} and is not available in the workflow.".format(
                    runname=runname, loc=loc), file=sys.stderr)
            elif not os.path.exists(os.path.join(loc, config['fast5_dir'])) or not os.listdir(os.path.join(loc, config['fast5_dir'])):
                print_("[WARNING] {runname} configured but with missing/empty reads directory.".format(
                    runname=runname), file=sys.stderr)
    if not config['nanopore']:
        print_("[WARNING] 'do_basecalling' set to True but 'nanopore' set to False. This will not end well.")
# check for sequences
else:
    for tag in config['runs']:
        sequences = os.path.join('sequences', tag+'.fastq.gz')

        if not os.path.exists(sequences):

            if config['merge_paired_end'] == True:
                if 'fwdReads' not in config['runs'][tag] or 'fwdReads' not in config['runs'][tag]:
                    print_(f"[WARNING] merge_paired_end set to True but forward and/or reverse reads files not provided for {tag} with keyword `fwdReads` and `rvsReads`")
                fwd = os.path.join('sequences', 'paired', config['runs'][tag]['fwdReads'])
                rvs = os.path.join('sequences', 'paired', config['runs'][tag]['rvsReads'])
                if not all((os.path.exists(fwd), os.path.exists(rvs))):
                    print_(f"[WARNING] merge_paired_end set to True but forward and/or reverse reads files provided for {tag}, {fwd}, {rvs} do not exist")

            elif 'runname' not in config['runs'][tag]:
                print_(f"[WARNING] `do_basecalling` set to False and runname director(y/ies) not set for tag `{tag}`, but sequences file `{sequences}` not found", file=sys.stderr)

            else:
                for runname in config['runs'][tag]['runname']:
                    batch = os.path.join('sequences', 'batches', runname)
                    if config['porechop']:
                        porechopBatch = batch + '_porechop'
                        if not os.path.exists(batch) and not os.path.exists(porechopBatch):
                            print_(f"[WARNING] `do_basecalling` set to False, sequences file `{sequences}` not found, basecalled sequence batch folder `{batch}` not found, and porechop sequence batch folder `{porechopBatch}` not found.", file=sys.stderr)
                    elif not os.path.exists(batch):
                        print_(f"[WARNING] `do_basecalling` set to False, `porechop` set to False, sequences file `{sequences}` not found, and basecalled sequence batch folder `{batch}` not found", file=sys.stderr)

# check reference sequences
for option in ['do_AA_analysis', 'auto_detect_longest_ORF']:
    if option not in config:
        raise RuntimeError(f"[ERROR] required boolean option `{option}` not present in config file.")
refSeqErrors = []
refSeqFastaFiles = []   # list files for all tags then check so files not being checked multiple times
for tag in config['runs']:
    if 'reference' not in config['runs'][tag]:
        refSeqErrors.append(f"[ERROR] No reference file provided for tag `{tag}")
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
    if config['do_AA_analysis'] == False:
        if 'AA_muts_of_interest' in config['runs'][tag]:
            print_(f'[WARNING] AA_muts_of_interest provided for {tag}, but `do_AA_analysis` set to False. AA muts of interest will not be evaluated.', file=sys.stderr)

for refFasta, alnFasta in refSeqFastaFiles:
    referenceSeqs = list(SeqIO.parse(refFasta, 'fasta'))
    if len(referenceSeqs) < 2:
        refSeqErrors.append(f"[ERROR] Reference file {refFasta} must be fasta formatted and contain more than 2 sequences. See example_working_directory/ref/exampleReferences.fasta for more information.")
        continue
    elif len(referenceSeqs) >3:
        print_(f"[WARNING] Reference file {refFasta} contains more than the maximum usable number of three sequences. Only the first three sequences will be used for alignment and analysis.")
    alignmentSeq, nucleotideSeq = referenceSeqs[0], referenceSeqs[1]

    if config['do_AA_analysis'] == True:
        if config['auto_detect_longest_ORF'] == True:
            try:
                longestORF = max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)', str(nucleotideSeq.seq.upper())), key = len) # taken from https://stackoverflow.com/questions/31757876/python-find-longest-orf-in-dna-sequence#31758161
                proteinSeq = nucleotideSeq
                proteinSeq.id = proteinSeq.id + '_prot'
                proteinSeq.seq = longestORF
            except ValueError:
                if len(referenceSeqs) < 3:
                    refSeqErrors.append(f"[ERROR] Could not auto detect an ORF and no protein sequence provided for {refFasta}.")
                    continue
                print_(f"[WARNING] `auto_detect_longest_ORF` set to True, but no standard ORF found in reference sequence `{refFasta}`. Will use provided protein sequence.", file=sys.stderr)
        else:
            if len(referenceSeqs) < 3:
                refSeqErrors.append(f"[ERROR] `do_AA_analysis` set to True but no protein sequence provided for {refFasta}")
            proteinSeq = referenceSeqs[2]
        if nucleotideSeq.seq.upper().find(proteinSeq.seq.upper()) == -1:
            refSeqErrors.append(f"[ERROR] Protein reference sequence `{proteinSeq.id}` of reference file `{refFasta}` is not a subsequence of nucleotide sequence {nucleotideSeq.id}")
        if len(proteinSeq.seq)%3 != 0:
            refSeqErrors.append(f"[ERROR] Protein reference sequence `{proteinSeq.id}` of reference file `{refFasta}` length not a multiple of 3, and therefore cannot be used as ORF")
        for i, nt in enumerate(str(proteinSeq.seq).upper()):
            if nt not in list("ATGCN"):
                refSeqErrors.append(f"[ERROR] Character {nt} at position {i} in reference sequence `{proteinSeq.id}` of reference file `{refFasta}` is not a canonical nucleotide")
    elif len(referenceSeqs) > 2:
        print_(f"[WARNING] `do_AA_analysis` set to False, but more than 2 sequences provided in reference file {refFasta}. AA analysis will not be performed.")

    if alignmentSeq.seq.upper().find(nucleotideSeq.seq.upper()) == -1:
        refSeqErrors.append(f"[ERROR] Nucleotide (second) sequence `{nucleotideSeq.id}` is not a subsequence of alignment (first) `{alignmentSeq.id}` sequence in reference file `{refFasta}`")
    for i, nt in enumerate(str(nucleotideSeq.seq).upper()):
        if nt not in list("ATGCN"):
            refSeqErrors.append(f"[ERROR] Character {nt} at position {i} in reference sequence `{nucleotideSeq.id}` of reference file `{refFasta}` is not a canonical nucleotide")
    for i, nt in enumerate(str(nucleotideSeq.seq).upper()):
        if nt not in list("ATGCN"):
            refSeqErrors.append(f"[ERROR] Character {nt} at position {i} in reference sequence `{alignmentSeq.id}` of reference file `{refFasta}` is not a canonical nucleotide")

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
        print_(f'Alignment reference .fasta file not found or is different from original reference .fasta file. Generating {alnFasta} from {refFasta}.', file=sys.stderr)
        with open(alnFasta, 'w') as fastaOut:
            first_record = next(SeqIO.parse(refFasta, 'fasta'))
            fastaOut.write(f'>{first_record.id}\n{first_record.seq.upper()}\n')

if len(refSeqErrors) > 0:
    for err in refSeqErrors:
        print_(err, file=sys.stderr)
    raise RuntimeError("Errors in reference sequences found. See above.")

# UMI checks
if 'UMI_consensus' not in config:
    raise RuntimeError("[ERROR] Required boolean option `UMI_consensus` not present in config file.")
if config['UMI_consensus'] == True:
    consensusCopyDict = {}      # dictionary that is used to reuse consensus sequences from a different tag if they are generated using the same files. Keys are tags, and values are tags whose consensus sequences will be used for downstream files for the key tag
    consensusFilesDict = {}     # nested dictionary that keeps track of the runnames and reference sequence used for each tag. keys, subkeys, and subsubkeys are the alignment sequence, a list of the runnames, and a list of the UMIs respectively, and values are the first tag encountered that uses those runnames and alignment sequences
    for tag in config['runs']:
        if 'UMI_contexts' not in config['runs'][tag]:
            print_(f"[WARNING] `UMI_contexts` set to True, but tag `{tag}` does not contain the required key `UMI_contexts`. UMI consensus steps will fail.")
        elif 'UMI_contexts' in config['runs'][tag]:
            refFasta = config['runs'][tag]['reference']
            alignmentSeq = list(SeqIO.parse(refFasta, 'fasta'))[0]
            for i, context in enumerate(config['runs'][tag]['UMI_contexts']):
                if str(alignmentSeq.seq).upper().find(context.upper()) == -1:
                    print_(f"[WARNING] UMI context {i+1} for tag `{tag}`, `{context}`, not found in reference `{alignmentSeq.id}` in fasta `{refFasta}`. UMI consensus steps will fail.")
            runnames = config['runs'][tag]['runname']
            runnames.sort()
            runnames = tuple(runnames)
            UMIcontexts = config['runs'][tag]['UMI_contexts']
            UMIcontexts.sort()
            UMIcontexts = tuple(UMIcontexts)
            if str(alignmentSeq.seq).upper() not in consensusFilesDict:
                consensusFilesDict[str(alignmentSeq.seq).upper()] = {runnames:{UMIcontexts:tag}}
                consensusCopyDict[tag] = tag
            else:
                if runnames not in consensusFilesDict[str(alignmentSeq.seq).upper()]:
                    consensusFilesDict[str(alignmentSeq.seq).upper()][runnames] = {UMIcontexts:tag}
                    consensusCopyDict[tag] = tag
                else:
                    if UMIcontexts not in consensusFilesDict[str(alignmentSeq.seq).upper()][runnames]:
                        consensusFilesDict[str(alignmentSeq.seq).upper()][runnames][UMIcontexts] = tag
                        consensusCopyDict[tag] = tag
                    else:
                        consensusTag = consensusFilesDict[str(alignmentSeq.seq).upper()][runnames][UMIcontexts]
                        print_(f"[NOTICE] Runname(s), alignment sequence, and UMI context combination used more than once. Using the consensus .fasta file of tag `{consensusTag}` for tag `{tag}` to reduce computation time and storage requirements")
                        consensusCopyDict[tag] = consensusTag
    config['consensusCopyDict'] = consensusCopyDict

# Demultiplexing checks
if 'demux' not in config:
    raise RuntimeError("[ERROR] Required boolean option `demux` not present in config file.")
if config['demux'] == True:
    for tag in config['runs']:
        if 'barcodeInfo' not in config['runs'][tag]:
            print_(f"[WARNING] `demux` set to True, but tag {tag} does not contain key `barcodeInfo`. Demultiplexing will fail.")
        if len(config['runs'][tag]['barcodeInfo']) == 0:
            print_(f"[WARNING] `demux` set to True, but tag {tag} `barcodeInfo` does not contain any barcode types. Demultiplexing will fail.")
        refFasta = config['runs'][tag]['reference']
        alignmentSeq = list(SeqIO.parse(refFasta, 'fasta'))[0]
        for barcodeType in config['runs'][tag]['barcodeInfo']:
            for requiredKey in ['context', 'fasta', 'reverseComplement']:
                if requiredKey not in config['runs'][tag]['barcodeInfo'][barcodeType]:
                    print_(f"[WARNING] `demux` set to True, but tag `{tag}` barcode type `{barcodeType}` does not contain the required key `{requiredKey}`.")
            if str(alignmentSeq.seq).upper().find(config['runs'][tag]['barcodeInfo'][barcodeType]['context'].upper()) == -1:
                print_(f"[WARNING] Barcode type `{barcodeType}` context `{config['runs'][tag]['barcodeInfo'][barcodeType]['context']}` not found in reference `{alignmentSeq.id}` in fasta `{refFasta}`")
            bcFasta = os.path.join(config['references_directory'], config['runs'][tag]['barcodeInfo'][barcodeType]['fasta'])
            config['runs'][tag]['barcodeInfo'][barcodeType]['fasta'] = bcFasta
            if os.path.isfile(bcFasta):
                if len(list(SeqIO.parse(bcFasta, 'fasta'))) == 0:
                    print_(f"[WARNING] Barcode fasta file `{bcFasta}` empty or not fasta format")
                if any(['_' in bc.id for bc in list(SeqIO.parse(bcFasta, 'fasta'))]):
                    print_(f"[WARNING] Sequence ID(s) in barcode fasta file `{bcFasta}` contain underscore(s), which may disrupt the pipeline. Please remove all underscores in sequence IDs.", file=sys.stderr)
                if type(config['runs'][tag]['barcodeInfo'][barcodeType]['reverseComplement'])!=bool:
                    print_(f"[WARNING] Barcode type {barcodeType} reverseComplement not bool (True or False)")
            if 'noSplit' in config['runs'][tag]['barcodeInfo'][barcodeType]:
                if config['runs'][tag]['barcodeInfo'][barcodeType]['noSplit'] == True:
                    for group in config['runs'][tag]['barcodeGroups']:
                        for bcType in config['runs'][tag]['barcodeGroups'][group]:
                            if bcType == barcodeType:
                                print_(f"[WARNING] `noSplit` set to True for barcode type `{barcodeType}` in run tag `{tag}`, but is used for naming in barcode group `{group}`. Demultiplexing will fail.")
        for group in config['runs'][tag]['barcodeGroups']:
            for bcType in config['runs'][tag]['barcodeGroups'][group]:
                if bcType not in config['runs'][tag]['barcodeInfo']:
                    print_(f"[WARNING] At least one barcode type in barcode group `{group}` for run tag `{tag}` is not defined in 'barcodeInfo'. Demultiplexing will fail.")
                if config['runs'][tag]['barcodeGroups'][group][bcType] not in [seq.id for seq in list(SeqIO.parse(bcFasta, 'fasta'))]:
                    print_(f"[WARNING] Barcode type `{bcType}` in barcode group `{group}` for run tag `{tag}` is not present in the barcode fasta file `{config['runs'][tag]['barcodeInfo'][bcType]['fasta']}` set for this tag")
            
elif config['demux'] == False:
    for tag in config['runs']:
        if ('barcodeInfo' or 'barcodeGroups') in config['runs'][tag]:
            print_(f"[WARNING] `barcodeInfo` or `barcodeGroups` provided for run tag `{tag}` but `demux` set to False. Demultiplexing will not be performed.", file=sys.stderr)

# check that tags and barcodeGroup names don't contain underscores
for tag in config['runs']:
    if '_' in tag:
        print_(f"[WARNING] Run tag `{tag}` contains underscore(s), which will disrupt the pipeline. Please remove all underscores in run tag names.", file=sys.stderr)
    if 'barcodeGroups' in config['runs'][tag]:
        for bcGroup in config['runs'][tag]['barcodeGroups']:
            if '_' in bcGroup:
                print_(f"[WARNING] Barcode group `{bcGroup}` for run tag `{tag}` contains underscore(s), which will disrupt the pipeline. Please remove all underscores in barcode group names.", file=sys.stderr)

# add timepoints files to config dictionary in the format {'timepoints':{tag:timepointCSVfile}}. Timepoint CSV files are not used more than once
#   and will instead be assigned to the first tag that uses that file. Also checks for the following:
#       - csv file exists
#       - Referenced tags and barcodeGroups are defined in config file
#       - reference fasta files for sample/barcodeGroup combinations are the same in each row
#       - at least two timepoints are given
#       - a row only uses tags from the same sequencing run or, if it uses different sequencing runs, that a 'background'
#           barcode group is provided in the config file. This is important because background subtraction is necessary
#           for accurate rate calculations, and sequencing error can of course differ from run to run.
timepointErrors = []
for tag in config['runs']:
    if 'timepoints' in config['runs'][tag]:
        CSVpath = os.path.join(config['references_directory'], config['runs'][tag]['timepoints'])
        if 'timepoints' not in config:
            config['timepoints'] = {}
        if CSVpath not in config['timepoints'].values():
            config['timepoints'][tag] = CSVpath
        if os.path.isfile(CSVpath):
            timepointsCSV = pd.read_csv(CSVpath, index_col=0, header=1)
            topRow = [x for x in pd.read_csv(CSVpath)).columns if 'Unnamed: ' not in x]
                if len(topRow) > 1:
                    print_(f"[NOTICE] More than one cell is filled in the top row of timepoint CSV file {str(snakemake.input.timepoints)}. Only the first cell in this row will be used for labeling outputs of mutation rate plots.", file=sys.stderr)
                elif len(topRow) == 0: 
                    print_(f"[NOTICE] No time unit provided in top row of timepoint CSV file {str(snakemake.input.timepoints)}. Default 'generations' will be used.", file=sys.stderr)
            if len(timepointsCSV.columns) <= 1:
                print_(f"[WARNING] Timepoints .CSV file for run tag `{tag}`, `{CSVpath}` does not have at least two timepoints. Timepoint-based snakemake rules will fail.", file=sys.stderr)
            else:
                rowIndex = 2    # start at 2 because first two rows are ignored with pd.read_csv call
                for _, row in timepointsCSV.iterrows():
                    rowIndex += 1
                    firstTP = timepointsCSV.columns[0]
                    firstTag = row[firstTP].split('_')[0]
                    if firstTag in config['runs']:
                        firstTagRefFasta = config['runs'][firstTag]['reference']
                        firstTagRefSeq = str(list(SeqIO.parse(firstTagRefFasta, 'fasta'))[1].seq).upper()
                        firstTagRunname = config['runs'][firstTag]['runname']
                    else:
                        timepointErrors.append(f"[ERROR] Tag referenced in row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{firstTag}` is not defined in config file. Check timepoints csv file and config file for errors.")
                    for tp in timepointsCSV.columns[1:]:
                        tag = row[tp].split('_')[0]
                        if tag in config['runs']:
                            if firstTag in config['runs']:
                                tagRefSeq = str(list(SeqIO.parse(config['runs'][tag]['reference'], 'fasta'))[1].seq).upper()
                                if tagRefSeq != firstTagRefSeq:
                                    print_(f"[WARNING] In row {rowIndex} of timepoints .CSV file `{CSVpath}`, samples `{row[firstTP]}` and `{row[tp]}` use different reference sequences. Analysis may be unreliable.", file=sys.stderr)
                                tagRunname = config['runs'][tag]['runname']
                                if tagRunname != firstTagRunname and 'background' not in config:
                                    print_(f"[WARNING] In row {rowIndex} of timepoints .CSV file `{CSVpath}`, samples `{row[firstTP]}` and `{row[tp]}` use different runnames, but a background barcodeGroup is not provided. Analysis may be unreliable.", file=sys.stderr)
                        else:
                            timepointErrors.append(f"[ERROR] Tag referenced in row {rowIndex} of timepoints .CSV file `{CSVpath}`, `{tag}` is not defined in config file. Check timepoints csv file and config file for errors.")
                    
        else:
            print_(f"[WARNING] Timepoints .CSV file for run tag `{tag}`, `{CSVpath}` does not exist.", file=sys.stderr)
if len(timepointErrors) > 0:
    for err in timepointErrors:
        print_(err, file=sys.stderr)
    raise RuntimeError("Errors in timepoint CSV files found. See above.")
        

# # include modules
include : "rules/storage.smk"
include : "rules/clean.smk"
include : "rules/report.smk"
include : "rules/pipeline.smk"

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

# ---------------------------------------------------------------------------------
# Copyright (c) 2018-2020, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann, modified by Gordon Rix
# ---------------------------------------------------------------------------------

def targets_input(wildcards):
    out = []
    out.append('mutation-stats.csv')
    if config['demux'] == True:
        out.append('demux-stats.csv')
    out.extend(expand('plots/{tag}_{AAorNT}-mutation-distributions.html', tag=config['runs'], AAorNT=['AA','NT'] if config['do_AA_analysis'] else ['NT']))
    out.extend(expand('plots/{tag}_mutation-spectra.html', tag=config['runs']))
    out.extend(expand('plots/{tag}_{AAorNT}-mutations-frequencies.html', tag=config['runs'], AAorNT=['AA','NT'] if config['do_AA_analysis'] else ['NT']))
    if config['nanoplot'] == True:
        out.extend(expand('plots/nanoplot/{tag}_fastq_NanoStats.txt', tag=config['runs']))
        out.extend(expand('plots/nanoplot/{tag}_alignment_NanoStats.txt', tag=config['runs']))
    if config['UMI_consensus'] == True:
        out.extend(expand('plots/{tag}_UMIgroup-distribution.html', tag=list(set( [config['consensusCopyDict'][str(t)] for t in config['runs']] )) ))
        if config['nanoplot'] == True:
            out.extend(expand('plots/nanoplot/{tag}_alignment_preConsensus_NanoStats.txt', tag=list(set( [config['consensusCopyDict'][str(t)] for t in config['runs']] ))))
        out.append('sequences/UMI/UMI-extract-summary.csv')
    if ('dms_view_chain' and 'dms_view_chain_numbering_difference') in config and config['do_AA_analysis'] == True:
        out.append('dms-view-table.csv')
    out.extend(expand('plots/{tag}_pipeline-throughput.html', tag=config['runs']))
    if 'timepoints' in config:
        out.extend(expand('plots/{tag}_mutation-rates.html', tag=config['timepoints']))
        # out.extend(expand('plots/{tag}_mutation-rate-spectrum.html', tag=config['timepoints']))

    if config['diversity_plot_all']:
        if config['demux']:
            out.extend( expand('plots/.{tag}_allDiversityPlots.done', tag=config['runs']))
        else:
            out.extend( expand('mutation_data/{tag}_all_{dataType}', tag=config['runs'], dataType =['diversity-graph.gexf', 'hamming-distance-distribution.csv']) )
            out.extend( expand('plots/{tag}_all_{plotType}', tag=config['runs'], plotType =['diversity-graph.html', 'hamming-distance-distribution.html']) )
    elif config.get('diversity_plot_subset', False) not in ['',False]:
        out.extend( expand('mutation_data/{tag_barcodes}_{dataType}', tag_barcodes=config['diversity_plot_subset'].split(', '), dataType=['diversity-graph.gexf', 'hamming-distance-distribution.csv']) )
        out.extend( expand('plots/{tag_barcodes}_{plotType}', tag_barcodes=config['diversity_plot_subset'].split(', '), plotType=['hamming-distance-distribution.html', 'diversity-graph.html']) )

    return out

rule targets:
    input:
        targets_input