# Maple: Mutation Analysis for Parallel Evolution

Maple is a snakemake-based workflow for analysis of mutation-rich next generation sequencing data
from highly parallelized targeted evolution experiments, with a special focus on Oxford Nanopore
Technology (ONT) data. It provides consensus sequence generation using unique molecular
identifiers (UMIs), easy-to-use and versatile demultiplexing, and a suite of analyses and plots.

## Getting started

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
