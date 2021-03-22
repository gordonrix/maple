# Maple: Mutation Analysis for Parallel Evolution

Maple is a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline for analysis of
mutation-rich next generation sequencing data from highly parallelized targeted evolution experiments, with a
special focus on Oxford Nanopore Technology (ONT) data. It provides consensus sequence generation using
unique molecular identifiers (UMIs), easy-to-use and versatile demultiplexing, and a suite of analyses and plots.

Analysis is primarily performed by custom python scripts, but several external tools are integrated
and critical:
 - [Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis)
 - [minimap2](https://doi.org/10.1093/bioinformatics/bty191)
 - [Samtools](http://www.htslib.org/)
 - [NGmerge](https://github.com/harvardinformatics/NGmerge)

Additionally, many concepts and code are borrowed from the (very sophisticated) snakemake pipeline [Nanopype](https://nanopype.readthedocs.io/en/latest/)

## Setup

Maple requires conda, which can be installed by following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
Miniconda is lighter weight and provides all that is needed by Maple. If basecalling ONT data is needed,
then GPU hardware and [cuda](https://docs.nvidia.com/cuda/) are also required. To accomodate use on
compute clusters, the only steps that require root privileges are installation of conda and cuda,
which should already be available on most clusters.

Clone the repository:

    git clone https://github.com/gordonrix/maple.git
    cd maple


Create a conda environment that enables usage of mamba, which is better capable of creating environments
with snakemake:

    conda create mamba -n mambaEnv -c conda-forge


With root privileges, mamba may instead be installed in the base environment to reduce package redundancy:

    conda install mamba -n base -c conda-forge


Next, use mamba to create the environment that we need. If using mamba from the base environment, the default
location is probably fine, but if using mamba within its own environment, the path to the environment should
be specified. In this example, the default path prefix ~/.conda/envs is used, but please ensure this is the correct location:

    conda activate mambaEnv
    mamba create --prefix ~/.conda/envs/maple -c bioconda -c conda-forge --file requirements.txt -y
    conda deactivate
    conda activate maple


The various tools that maple uses that are not simple python scripts must now be installed. This is accomplished
easily through snakemake, using the 'install' snakefile, derived from a similar setup in [Nanopype](https://nanopype.readthedocs.io/en/latest/).
Running the following command, with the --prefix flag modified to point to the location of the environment
we created in the previous step, will carry out all the steps to provide the environment with program binaries
needed for execution of the pipeline:

    snakemake --snakefile rules/install.smk --directory ~/.conda/envs/maple -j 4 all


## Usage

Snakemake takes as input a file name corresponding to a file that it is capable of generating via some sequence of steps, or 'rules', working backwards to determine
which rules must be run, then running each of these rules to generate the user-specified file, along with all the files that were generated while carrying out
the necessary steps.

Maple relies upon a configuration file to determine what files can be generated. A portion of this config file is devoted to global settings that will
be applied to all datasets analyzed using this config file, while another portion is used to designate information specific to each dataset being analyzed,
organized by 'tags' which will be applied to all filenames. Detailed information on these settings can be found in example_working_directory/config.yaml.

Following installation, the steps required for basic usage of the pipeline are as follows:
1. Create a directory for the dataset(s) that you wish to analyze by copying the example_working_directory to a location of your choice,
rename it as desired, and navigate to this new directory.
2. Modify the config.yaml file within the new directory as appropriate. Most of these settings can remain unchanged, but the settings for individual
tag(s) should be changed, such as the sequencing data and reference sequences file locations and barcode information (if demultiplexing is required)
3. Modify the reference sequence and barcode .fasta files (located in the 'ref' directory) as appropriate. Use the exampleReferences.fasta file
for assistance with determining the appropriate reference sequences.
4. Activate the maple conda environment that you created during installation if it's not already active, and run snakemake by requesting a specific file,
and designating a number of threads. In most cases, at least 4 threads should be used. Take care to run the pipeline only when in the working
directory (e.g. example_working_directory), otherwise the --directory flag must be used to specify a directory that contains the appropriate
files in the correct relative locations (i.e. config.yaml, ref/*, etc.). The path to the maple snakefile must also be modified as appropriate:

    conda activate maple
    snakemake --snakefile PATH/TO/maple/Snakefile -j 4 example_mutation-stats.csv

In place of a specific file name, 'targets' can be used to invoke a rule that automatically carries out most of the analysis that maple can do
for each of the designated tags:

    snakemake --snakefile PATH/TO/maple/Snakefile -j 4 targets

Likewise, if you'd like to restart your analysis without cluttering your working directory with additional tags, or if you just want to package up the key analysis files
for transfer or storage, the 'clean' rule can be called. This will move or copy all the small files generated during analyses to a timestamped directory
and removes large files such as alignment files, without modifying large important files such as .fast5 files and unmerged .fastq files.

    snakemake --snakefile PATH/TO/maple/Snakefile -j 4 clean