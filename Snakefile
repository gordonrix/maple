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


# verify given references
for tag in config['runs']:
    if config['runs'][tag]['reference']:
        refName = config['runs'][tag]['reference']
        if not os.path.isfile(refName):
            print_(f'[WARNING] Reference .fasta file for {refName} (given path: {config[refName]}) not found.', file=sys.stderr)
    else:
        print_(f"[WARNING] No reference .fasta file provided for {tag} in config.yaml. Alignment and downstream tools will not work.", file=sys.stderr)


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
# check for sequences
else:
    for tag in config['runs']:
        sequences = os.path.join('sequences', tag+'.fastq.gz')

        if config['merge_paired_end'] == True:
            if not os.path.exists(sequences):
                if 'fwdReads' not in config['runs'][tag] or 'fwdReads' not in config['runs'][tag]:
                    print_(f"[WARNING] merge_paired_end set to True but forward and/or reverse reads files not provided for {tag} with keyword `fwdReads` and `rvsReads`")
                fwd = os.path.join('sequences', 'paired', config['runs'][tag]['fwdReads'])
                rvs = os.path.join('sequences', 'paired', config['runs'][tag]['rvsReads'])
                if not all((os.path.exists(fwd), os.path.exists(rvs))):
                    print_(f"[WARNING] merge_paired_end set to True but forward and/or reverse reads files provided for {tag}, {fwd}, {rvs} do not exist")

        elif 'runname' not in config['runs'][tag]:
            if not os.path.exists(sequences):
                print_(f"[WARNING] `do_basecalling` set to False and runname director(y/ies) not set for {config['runs'][tag]}, but sequences file `{sequences}` not found", file=sys.stderr)
        else:
            if not os.path.exists(sequences):
                for runname in config['runs'][tag]['runname']:
                    batch = os.path.join('sequences', 'batches', runname)
                    if not os.path.exists(batch):
                        print_(f"[WARNING] `do_basecalling` set to False, sequences file `{sequences}` not found, and basecalled sequence batch folder `{batch}` not found", file=sys.stderr)

# check reference sequences
for tag in config['runs']:
    refFasta = config['runs'][tag]['reference']
    referenceSeqs = list(SeqIO.parse(refFasta, 'fasta'))
    if len(referenceSeqs) == 3:
        alignmentSeq, nucleotideSeq, proteinSeq = referenceSeqs
        assert nucleotideSeq.seq.upper().find(proteinSeq.seq.upper()) != -1, f"Protein (third) sequence is not a subsequence of nucleotide (second) sequence in reference file `{refFasta}`"
        for i, nt in enumerate(str(proteinSeq.seq).upper()):
            assert nt in list("ATGCN"), f"Position {i} in reference sequence `{proteinSeq.id}` is not a canonical nucleotide"
    elif len(referenceSeqs) == 2:
        alignmentSeq, nucleotideSeq = referenceSeqs
    else:
        raise RuntimeError("Reference fasta file must be fasta formatted and contain 2 or 3 sequences. See example_working_directory/ref/exampleReferences.fasta for more information.")
    assert nucleotideSeq.seq.upper().find(proteinSeq.seq.upper()) != -1, f"Nucleotide (second) sequence is not a subsequence of alignment (first) sequence in reference file `{refFasta}`"
    for i, nt in enumerate(str(nucleotideSeq.seq).upper()):
        assert nt in list("ATGCN"), f"Position {i} in reference sequence `{nucleotideSeq.id}` is not a canonical nucleotide"
    for i, nt in enumerate(str(nucleotideSeq.seq).upper()):
        assert nt in list("ATGCN"), f"Position {i} in reference sequence `{alignmentSeq.id}` is not a canonical nucleotide"

# ADD CHECKS FOR UMIS.

# ADD A CHECK FOR BARCODE INFO. If 'barcodeInfo' or 'barcodeGroups' is present in {tag}, both must be present and all barcode types in barcode groups
# must be defined in 'barcodeInfo' dict. Probably did this in the beginning of the demux script, so can move that here.
if config['demux']:
    for tag in config['runs']:
        if 'barcodeInfo' not in config['runs'][tag]:
            raise RuntimeError(f"[ERROR] `demux` set to True, but tag {tag} does not contain key `barcodeInfo`.")
        if len(config['runs'][tag]['barcodeInfo']) == 0:
            raise RuntimeError(f"[ERROR] `demux` set to True, but tag {tag} `barcodeInfo` does not contain any barcode types.")
        refFasta = config['runs'][tag]['reference']
        alignmentSeq = list(SeqIO.parse(refFasta, 'fasta'))[0]
        for barcodeType in config['runs'][tag]['barcodeInfo']:
            for requiredKey in ['context', 'fasta', 'reverseComplement']:
                if requiredKey not in config['runs'][tag]['barcodeInfo'][barcodeType]:
                    raise RuntimeError(f"[ERROR] `demux` set to True, but tag {tag} barcode type `{barcodeType}` does not contain the required key `{requiredKey}`.")
            assert str(alignmentSeq.seq).upper().find(config['runs'][tag]['barcodeInfo'][barcodeType]['context'].upper()) != -1, f"Barcode type `{barcodeType}` context `{config['runs'][tag]['barcodeInfo'][barcodeType]['context']}` not found in reference `{alignmentSeq.id}` in fasta `{refFasta}`"
            bcFasta = config['runs'][tag]['barcodeInfo'][barcodeType]['fasta']
            if os.path.isfile(bcFasta):
                assert len(list(SeqIO.parse(bcFasta, 'fasta'))) != 0, f"Barcode fasta file `{bcFasta}` empty or not fasta format"
                assert type(config['runs'][tag]['barcodeInfo'][barcodeType]['reverseComplement'])==bool, f"Barcode type {barcodeType} reverseComplement not bool (True or False)"
else:
    for tag in config['runs']:
        if ('barcodeInfo' or 'barcodeGroups' in config['runs'][tag]) and config['demux']==False:
            print_(f"[WARNING] `barcodeInfo` or `barcodeGroups` provided for run tag `{tag}` but `demux` set to False. Demultiplexing will not be performed.", file=sys.stderr)
        

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
    out.extend(expand('{tag}_mutation-stats.csv', tag=config['runs']))
    out.extend(expand('plots/{tag}_{AAorNT}-mutation-distributions.html', tag=config['runs'], AAorNT=['AA','NT'] if config['do_AA_analysis'] else ['NT']))
    out.extend(expand('plots/{tag}_mutation-spectra.html', tag=config['runs']))
    out.extend(expand('plots/{tag}_{AAorNT}-mutations-frequencies.html', tag=config['runs'], AAorNT=['AA','NT'] if config['do_AA_analysis'] else ['NT']))
    if config['UMI_consensus']:
        out.extend(expand('plots/{tag}_UMIgroup-distribution.html', tag=config['runs']))
    return out

rule targets:
    input:
        targets_input