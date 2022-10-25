# Maple: Mutation Analysis for Parallel Evolution

Maple is a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline for analysis of
mutation-rich next generation sequencing data from highly parallelized targeted evolution experiments, with a
special focus on Oxford Nanopore Technology (ONT) data. It provides consensus sequence generation using both concatemer-
and unique molecular identifiers (UMI)-based consensus, easy-to-use and versatile demultiplexing, and a suite of analyses and plots.

Analysis is primarily performed by a mix of custom python scripts and several external tools:
 - [medaka](https://github.com/nanoporetech/medaka)
 - [minimap2](https://doi.org/10.1093/bioinformatics/bty191)
 - [Samtools](http://www.htslib.org/)
 - [NGmerge](https://github.com/harvardinformatics/NGmerge)
 - [NanoPlot](https://github.com/wdecoster/NanoPlot)
 - [C3POa](https://github.com/christopher-vollmers/C3POa)

Additionally, many concepts and code are borrowed from the snakemake pipeline [Nanopype](https://nanopype.readthedocs.io/en/latest/)

## Setup

Maple requires conda, which can be installed by following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
Miniconda is lighter weight and provides all that is needed by Maple. If basecalling ONT data is needed,
then GPU hardware and [cuda](https://docs.nvidia.com/cuda/) are also required. CPU basecalling is not available
by default because it is prohibitively slow. If CPU basecalling is desired, or if basecalling is not needed,
comment/uncomment the appropriate rules in rules/install.smk

Additionally, the following software packages are required:
 - git gcc g++ wget rsync
 - binutils autoconf make cmake
 - libgomp1
 - zlib1g-dev
 - bzip2 libbz2-dev
 - liblzma-dev libncurses5-dev
 - libcunit1 libhdf5-100 libidn11 libopenblas-base
 - libgssapi-krb5-2
 - libzstd-dev
 
These are likely already present in most production environments.

Clone the repository:

    git clone https://github.com/gordonrix/maple.git
    cd maple


Create a conda environment that enables usage of mamba, which is better capable of creating environments
with snakemake:

    conda create mamba -n mambaEnv -c conda-forge


Alternatively, with root privileges, mamba may instead be installed in the base environment to reduce package redundancy:

    conda install mamba -n base -c conda-forge


Next, use mamba to create the environment that we need. If using mamba within its own environment, the path to the environment should
be specified with the `--prefix In this example, the default path prefix ~/miniconda3/envs is used, but please ensure this is the correct location:

    conda activate mambaEnv
    mamba env create --prefix ~/miniconda3/envs/maple --file requirements.yaml
    conda deactivate
    conda activate maple


The various tools that maple uses that are not simple python scripts must now be installed. This is accomplished
easily through snakemake, using the 'install' snakefile, derived from a similar setup in [Nanopype](https://nanopype.readthedocs.io/en/latest/).
Running the following command, with the --prefix flag modified to point to the location of the environment
we created in the previous step, will carry out all the steps to provide the environment with program binaries
needed for execution of the pipeline:

    snakemake --snakefile rules/install.smk --directory ~/miniconda3/envs/maple -j 4 all

## Usage

Snakemake takes as input a file name corresponding to a file that it is capable of generating via some sequence of steps, or 'rules', working backwards to determine
which rules must be run, then running each of these rules to generate the user-specified file, along with all the files that were generated while carrying out
the necessary steps.

Maple relies upon a configuration file to determine what files can be generated. A portion of this config file is devoted to global settings that will
be applied to all datasets analyzed using this config file, while another portion is used to designate information specific to each dataset being analyzed,
organized by 'tags' which will be applied to all filenames. Detailed information on these settings can be found in example_working_directory/config.yaml.

# Running the example
Included in the github repository is an example_working_directory that will allow you to test that everything was properly installed, and will serve as a
jumping off point for future working directories. To test that everything is working, navigate to this working directory and perform a 'dry-run' (-n) for 
all files that are demanded by the `targets` directive to see that snakemake knows how to perform all the required steps.
Here we also include the '-q' or '--quiet' flag so that any warnings that might come up don't require much scrolling to get to:

    snakemake --snakefile ../Snakefile -j 4 targets -n -q

This should produce some gold text that shows all the different jobs that will be run. Running the same command without the extra flags will run all the rules:

    snakemake --snakefile ../Snakefile -j 4 targets

If this worked, you will see a stream of green text describing the steps that are being performed, followed by the text
    maple completed successfully

You can then explore the analysis files that have been created, or apply this same process to your own data.

# Running maple on your data

Following installation and running the example, the steps required for basic usage of the pipeline are as follows:
1. Create a directory for the dataset(s) that you wish to analyze by creating a new working directory with the name of your choice (preferably related to
a specific experiment) and copy the 'ref' directory and 'config.yaml' file from the example_working_directory to this new directory.
2. Modify the config.yaml file within the new directory as appropriate. Global settings (applied to the working directory using this config file) are above
the 'runs' key. Beneath the 'runs' key you must define the analysis for each individual run 'tag' that you wish to analyze. Ideally, a static directory that
contains all your nanopore sequencing data will be chosen and can be defined at the top of the config.yaml such that this setting need not be changed
when you copy the config.yaml to future working directories. It's not imperative to change most of the other global settings, though growing familiar
with them will be useful. Most of the settings for the 'tag'(s) you are analyzing should be changed, such as the sequencing data type/location and the
different required inputs for different types of analysis that should be performed (e.g. concatemer/UMI consensus, demultiplexing).
3. Modify the reference sequence and barcode .fasta files (located in the 'ref' directory) as appropriate. Use the exampleReferences.fasta file
for assistance with determining the appropriate reference sequences.
4. Activate the maple conda environment that you created during installation if it's not already active, and run snakemake by requesting a specific file,
and designating a number of threads. In most cases, at least 4 threads should be used. Take care to run the pipeline only when in the working
directory (e.g. example_working_directory), otherwise the --directory flag must be used to specify a directory that contains the appropriate
files in the correct relative locations (i.e. config.yaml, ref/*, etc.). The path to the maple snakefile must also be modified as appropriate.
Here I will ask maple to produce the mutation stats summary file for the TrpB tag, which is defined in the config file:

    conda activate maple
    snakemake --snakefile PATH/TO/maple/Snakefile -j 4 TrpB_mutation-stats.csv

Use of the '-n' flag is strongly recommended prior to running the full pipeline. This causes snakemake to do a 'dry-run' in which jobs are planned out, but
not executed. Because many checks are performed to identify any potential problems in how things were set up (e.g. checking that reference files
exist), this will better guarantee that the entire pipeline will run to completion prior to starting it. The '-q' flag should also be used so
that any warnings can be easily read.

In place of a specific file name, 'targets' can be used to invoke a rule that automatically carries out most of the analysis that maple can do
for each of the designated tags, as was done above:

    snakemake --snakefile PATH/TO/maple/Snakefile -j 4 targets

Likewise, if you'd like to restart your analysis without cluttering your working directory with additional tags, or if you just want to package up the key analysis files
for transfer or storage, the 'clean' rule can be called. This will move or copy all the small files generated during analyses to a timestamped directory
and removes large files such as alignment files, without modifying large important files such as .fast5 files and .fastq files. If analysis rules that
produce files in the 'sequences' directory need to be rerun, such as UMI rules or paired end merging, the outputs of those rules must be manually deleted or renamed
to enable rule re-run.

    snakemake --snakefile PATH/TO/maple/Snakefile -j 4 clean
