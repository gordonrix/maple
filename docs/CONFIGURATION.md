# Configuration Reference

This document provides a comprehensive reference for all configuration parameters in Maple. They are defined in the `config.yaml` file, which should be found at the top level of the Maple working directory.

## Required Configuration

<details>
<summary>metadata</summary>

Directory name for metadata files, located in working directory.

**Type:** String (directory name)  
**Default:** `metadata`  
**Example:** `metadata: metadata`

This directory should contain all files that define a Maple run except for the config file and sequencing data. For example, CSV files defining runs or barcode info and fasta files defining barcodes

</details>

<details>
<summary>runs</summary>

Path to CSV file defining experimental runs to be analyzed.

**Type:** String (file path)  
**Default:** None (must be specified)  
**Example:** `runs: tags.csv`

The tags.csv file must be located within the `metadata` directory and contains columns defining each experimental run:
- `tag`: A unique identifier for the run that will be used in all output filenames (**Note:** Tags may not contain underscores)
- `reference`: Reference FASTA filename in the metadata directory
- one of the following for sequence import:
    - `runname`: Directory or individual file name containing FASTQ files/sequences, must be within the config-defined `sequences_dir`
    - `bs_project_ID` and `sample_ID`: For pulling pre-demultiplexed data from Illumina's BaseSpace. Downloads the specified project ID and grabs the fastq.gz for the specified sample ID.
    - `fwdReads` and `rvsReads`: For paired-end read merging
- optional columns:
    - `barcode_info_csv`: CSV file defining barcode types and contexts
    - `partition_barcode_groups_csv`: CSV file defining names to assign to barcode groups for demultiplexing and partitioning sequences to distinct files
    - `label_barcode_groups_csv`: CSV file defining names to assign to barcode groups for demultiplexing and labeling sequences (in the 'BC' tag of a sequence in BAM file output)
    - `splint`: Splint sequence for RCA consensus generation
    - `UMI_contexts`: Context patterns for UMI extraction, in order of appearance in the reference, separated by `__` (double underscore)
    - `timepoint`: A CSV file used for time series analysis (see below)

</details>

<details>
<summary>sequences_dir</summary>

Path to directory containing raw sequencing data.

**Type:** String (directory path)  
**Default:** `data`  
**Example:** `sequences_dir: /path/to/sequences`

This directory will be searched for subdirectories matching the runname values from tags.csv. Can be an absolute path to a central sequence storage location that doesn't need to change between working directories.

</details>

<details>
<summary>fastq_dir (required)</summary>

Comma-separated list of folders within run directories to pull FASTQ files from.

**Type:** String (comma-separated)  
**Default:** `fastq_pass, fastq_fail`  
**Example:** `fastq_dir: fastq_pass`

Common configurations:
- `fastq_pass`: Only high-quality reads (Nanopore)
- `fastq_pass, fastq_fail`: All reads (Nanopore)
- `raw`: Custom directory name

</details>

## Consensus Generation

<details>
<summary>RCA Consensus Parameters</summary>

Parameters for rolling circle amplification consensus using C3POa.

**peak_finder_settings**
- **Type:** String
- **Default:** `'23,3,27,2'`
- **Description:** Settings used to identify splint alignment locations for splitting reads into subreads
- **Example:** `peak_finder_settings: '23,3,27,2'`

**RCA_batch_size**
- **Type:** Integer
- **Default:** `10000`
- **Description:** Number of sequences to batch together in each subprocess. Lower this number if RCA processing crashes
- **Example:** `RCA_batch_size: 5000`

**RCA_consensus_minimum**
- **Type:** Integer
- **Default:** `3`
- **Description:** Inclusive minimum number of complete subreads required to generate an RCA consensus read
- **Example:** `RCA_consensus_minimum: 3`

**RCA_consensus_maximum**
- **Type:** Integer
- **Default:** `20`
- **Description:** Inclusive maximum number of complete subreads used to generate an RCA consensus read. Reads with more subreads will not be used
- **Example:** `RCA_consensus_maximum: 15`

</details>

<details>
<summary>UMI Consensus Parameters</summary>

Parameters for UMI-based consensus generation using medaka and umicollapse.

**UMI_mismatches**
- **Type:** Integer
- **Default:** None
- **Description:** Maximum allowable number of mismatches that UMIs can contain and still be grouped together
- **Example:** `UMI_mismatches: 1`

**UMI_consensus_minimum**
- **Type:** Integer
- **Default:** None
- **Description:** Inclusive minimum number of subreads required to generate a UMI consensus read
- **Example:** `UMI_consensus_minimum: 5`

**UMI_consensus_maximum**
- **Type:** Integer
- **Default:** None
- **Description:** Inclusive maximum number of subreads used to generate a UMI consensus read. Groups with more subreads will be downsampled
- **Example:** `UMI_consensus_maximum: 20`

**UMI_medaka_batches**
- **Type:** Integer
- **Default:** None
- **Description:** Number of files to split BAM file into prior to running medaka. Increase if medaka throws memory errors
- **Example:** `UMI_medaka_batches: 10`

</details>

<details>
<summary>Medaka Parameters</summary>

Parameters for medaka consensus generation (Nanopore only).

**medaka_model**
- **Type:** String
- **Default:** `'r104_e81_sup_variant_g610'`
- **Description:** Model for medaka to use. Use `medaka smolecule --help` to see all available options
- **Example:** `medaka_model: 'r103_sup_variant_g507'`

**medaka_flags**
- **Type:** String
- **Default:** `'--quiet'`
- **Description:** Additional flags to add to medaka smolecule command. Threads, chunk-len, and model flags are already added
- **Example:** `medaka_flags: '--quiet --chunk-len 10000'`

</details>

## Alignment and Processing

<details>
<summary>Alignment Parameters</summary>

Parameters for sequence alignment using minimap2 and samtools.

**alignment_minimap2_flags**
- **Type:** String or Dictionary
- **Default:** `'-a -A2 -B4 -O4 -E2 --end-bonus=30 --secondary=no'`
- **Description:** Command line flags for minimap2 DNA alignment. Default options are optimized for targeted sequencing
- **Example (string):** `alignment_minimap2_flags: '-a -A2 -B4 -O4 -E2 --secondary=no'`
- **Example (dict):** 
  ```yaml
  alignment_minimap2_flags:
    tag1: '-a -A2 -B4 -O4 -E2 --secondary=no'
    tag2: '-a -A2 -B4 -O10 -E4 --secondary=no'
  ```

**alignment_samtools_flags**
- **Type:** String
- **Default:** `''`
- **Description:** Additional flags for samtools during alignment processing
- **Example:** `alignment_samtools_flags: '-F 4'`

</details>

<details>
<summary>Thread Configuration</summary>

Thread allocation for different pipeline steps.

**threads_medaka**
- **Type:** Integer
- **Default:** `2`
- **Description:** Threads per medaka processing batch
- **Example:** `threads_medaka: 4`

**threads_alignment**
- **Type:** Integer
- **Default:** `4`
- **Description:** Threads for alignment step. Recommended to provide maximum available threads, minimum 3
- **Example:** `threads_alignment: 8`

**threads_samtools**
- **Type:** Integer
- **Default:** `1`
- **Description:** Threads for samtools operations
- **Example:** `threads_samtools: 2`

**threads_demux**
- **Type:** Integer
- **Default:** `4`
- **Description:** Threads for demultiplexing operations
- **Example:** `threads_demux: 8`

</details>

## Demultiplexing

<details>
<summary>Demultiplexing Parameters</summary>

Parameters controlling barcode detection and sequence sorting.

**demux_screen_no_group**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Set to True if sequences not assigned to a named barcode group should be blocked from subsequent analysis
- **Example:** `demux_screen_no_group: False`

**demux_screen_failures**
- **Type:** Boolean
- **Default:** `False`
- **Description:** Set to True if sequences that fail barcode detection should be blocked from subsequent analysis
- **Example:** `demux_screen_failures: True`

**demux_threshold**
- **Type:** Float
- **Default:** `0.01`
- **Description:** Minimum proportion of total reads required for a demultiplexed file to be processed further
- **Example:** `demux_threshold: 0.05`

</details>

## Paired-End Read Processing

<details>
<summary>NGmerge Parameters</summary>

Parameters for paired-end read merging using NGmerge.

**merge_paired_end**
- **Type:** Boolean
- **Default:** `False`
- **Description:** Set to True if merging of paired-end reads is needed
- **Example:** `merge_paired_end: True`

**NGmerge_flags**
- **Type:** String
- **Default:** `'-m 10'`
- **Description:** Command line flags for NGmerge. `-m X` sets minimum allowable overlap to X
- **Example:** `NGmerge_flags: '-m 15 -u 41 -g'`

</details>

## Mutation Analysis

<details>
<summary>Mutation Analysis Parameters</summary>

Parameters controlling mutation detection and analysis.

**mutation_analysis_quality_score_minimum**
- **Type:** Integer
- **Default:** `5`
- **Description:** Minimum quality score needed for mutation to be counted. For amino acid analysis, all nucleotides in the codon must be above threshold
- **Example:** `mutation_analysis_quality_score_minimum: 10`

**sequence_length_threshold**
- **Type:** Float
- **Default:** `0.1`
- **Description:** Proportion of sequence length used as threshold for discarding aberrant-length sequences
- **Example:** `sequence_length_threshold: 0.2`

**analyze_seqs_with_indels**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Set to True if sequences containing insertions or deletions should be analyzed
- **Example:** `analyze_seqs_with_indels: False`

**mutations_frequencies_raw**
- **Type:** Boolean
- **Default:** `False`
- **Description:** If True, outputs mutation frequencies as raw counts instead of proportions
- **Example:** `mutations_frequencies_raw: True`

**mutations_frequencies_number_of_positions**
- **Type:** Integer
- **Default:** `20`
- **Description:** Number of mutations to include in most/least frequent mutations plots
- **Example:** `mutations_frequencies_number_of_positions: 30`

**mutations_frequencies_heatmap**
- **Type:** Boolean
- **Default:** `True`
- **Description:** If True, frequencies plot will be a heatmap; otherwise a stacked bar chart
- **Example:** `mutations_frequencies_heatmap: False`

</details>

<details>
<summary>Genotype Analysis Parameters</summary>

Parameters for genotype identification and analysis.

**highest_abundance_genotypes**
- **Type:** Integer
- **Default:** `10`
- **Description:** Number of most frequently appearing genotypes to find representative sequences for
- **Example:** `highest_abundance_genotypes: 20`

**genotype_ID_alignments**
- **Type:** Integer or String
- **Default:** `0`
- **Description:** Comma-separated list of genotype IDs to include in output, or 0 if not desired
- **Example:** `genotype_ID_alignments: '1,5,10'`

**unique_genotypes_count_threshold**
- **Type:** Integer
- **Default:** `5`
- **Description:** Minimum number of reads of a genotype for it to be included in unique genotypes count
- **Example:** `unique_genotypes_count_threshold: 10`

</details>

## Plotting and Visualization

<details>
<summary>Plot Configuration</summary>

Parameters controlling plot generation and appearance.

**colormap**
- **Type:** String
- **Default:** `'kbc_r'`
- **Description:** Colormap for plots. Options include 'kbc', 'fire', 'bgy', 'bgyw', 'bmy', 'gray', 'rainbow4', and their reverse versions with '_r'
- **Example:** `colormap: 'fire'`

**distribution_x_range**
- **Type:** Boolean or String
- **Default:** `False`
- **Description:** Comma-separated pair of values for x-axis range of distribution plots, or False for auto-range
- **Example:** `distribution_x_range: '0,100'`

**distribution_y_range**
- **Type:** Boolean or String
- **Default:** `False`
- **Description:** Comma-separated pair of values for y-axis range of distribution plots, or False for auto-range
- **Example:** `distribution_y_range: '0,0.5'`

**export_SVG**
- **Type:** Boolean or String
- **Default:** `False`
- **Description:** Controls SVG export. False=no export, True=export all, string=export plots containing string
- **Example:** `export_SVG: 'mutation'`

</details>

<details>
<summary>Plot Enable/Disable</summary>

Boolean flags to enable or disable specific plot types.

**plot_demux**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Generate demultiplexing statistics plots
- **Example:** `plot_demux: False`

**plot_mutation-distribution**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Generate mutation distribution plots
- **Example:** `plot_mutation-distribution: False`

**plot_mutation-frequencies**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Generate mutation frequency plots
- **Example:** `plot_mutation-frequencies: False`

**plot_hamming-distance-distribution**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Generate hamming distance distribution plots
- **Example:** `plot_hamming-distance-distribution: False`

**plot_mutation-distribution-violin**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Generate violin plots for mutation distributions
- **Example:** `plot_mutation-distribution-violin: False`

**plot_genotypes2D**
- **Type:** Boolean
- **Default:** `False`
- **Description:** Generate 2D genotype plots (also depends on other genotypes2D options)
- **Example:** `plot_genotypes2D: True`

</details>

## Advanced Features

<details>
<summary>Enrichment Analysis</summary>

Parameters for enrichment score calculation and filtering.

**enrichment_SE_filter**
- **Type:** Float
- **Default:** `0`
- **Description:** Proportion of standard errors to filter out (0-1). 0 or 1 disables filter
- **Example:** `enrichment_SE_filter: 0.1`

**enrichment_t0_filter**
- **Type:** Float
- **Default:** `0`
- **Description:** Proportion of timepoint 0 counts to filter out (0-1). 0 or 1 disables filter
- **Example:** `enrichment_t0_filter: 0.1`

**enrichment_score_filter**
- **Type:** Boolean or Float
- **Default:** `False`
- **Description:** Enrichment score threshold. Scores below this value are removed. False disables filter
- **Example:** `enrichment_score_filter: 0.5`

**enrichment_missing_replicates_filter**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Whether to filter out barcodes without enrichment scores for all replicates
- **Example:** `enrichment_missing_replicates_filter: False`

**enrichment_reference_bc**
- **Type:** String
- **Default:** `'all_barcodes'`
- **Description:** Barcode to use as reference for normalization. 'all_barcodes' normalizes to all barcodes
- **Example:** `enrichment_reference_bc: 'control'`

</details>

<details>
<summary>Genotypes2D Plotting</summary>

Parameters for 2D genotype visualization using dimensionality reduction.

**genotypes2D_plot_all**
- **Type:** Boolean
- **Default:** `False`
- **Description:** Generate 2D genotype plots for all samples individually
- **Example:** `genotypes2D_plot_all: True`

**genotypes2D_plot_groups**
- **Type:** Boolean
- **Default:** `False`
- **Description:** Generate 2D genotype plots for groups of samples
- **Example:** `genotypes2D_plot_groups: True`

**genotypes2D_plot_downsample**
- **Type:** Integer
- **Default:** `10000`
- **Description:** Maximum number of genotypes to include in 2D plots (downsampling threshold)
- **Example:** `genotypes2D_plot_downsample: 5000`

**genotypes2D_plot_AA**
- **Type:** Boolean
- **Default:** `True`
- **Description:** Use protein sequence for dimension reduction if available
- **Example:** `genotypes2D_plot_AA: False`

**genotypes2D_plot_point_size_col**
- **Type:** String
- **Default:** `count`
- **Description:** Genotypes column to use for point size (must be numerical)
- **Example:** `genotypes2D_plot_point_size_col: 'frequency'`

**genotypes2D_plot_point_size_range**
- **Type:** String
- **Default:** `'30, 60'`
- **Description:** Comma-separated pair of integers for minimum and maximum point sizes
- **Example:** `genotypes2D_plot_point_size_range: '20, 80'`

</details>

<details>
<summary>Hamming Distance Analysis</summary>

Parameters for hamming distance distribution analysis.

**hamming_distance_distribution_downsample**
- **Type:** Integer or Boolean
- **Default:** `1000`
- **Description:** Maximum number of genotypes to use for hamming distance calculation. False uses all genotypes
- **Example:** `hamming_distance_distribution_downsample: 500`

**hamming_distance_distribution_raw**
- **Type:** Boolean
- **Default:** `False`
- **Description:** If True, y-axis shows raw counts instead of proportions
- **Example:** `hamming_distance_distribution_raw: True`

</details>

<details>
<summary>Dashboard Configuration</summary>

Parameters for the interactive dashboard.

**dashboard_input**
- **Type:** String
- **Default:** `TrpB`
- **Description:** Tag, timepoint name, or tag_barcodeGroup to use for the dashboard
- **Example:** `dashboard_input: 'experiment1'`

**dashboard_port**
- **Type:** Integer
- **Default:** `3366`
- **Description:** Port number for running the dashboard
- **Example:** `dashboard_port: 8080`

</details>

<details>
<summary>Time Series Analysis</summary>

Parameters for time series and evolution experiments.

**timepoints**
- **Type:** String
- **Default:** `timepoints.csv`
- **Description:** CSV file providing tag and barcode combinations for experiment timepoints
- **Example:** `timepoints: my_timepoints.csv`

**uniques_only**
- **Type:** Boolean
- **Default:** `False`
- **Description:** If True, only uses unique mutations to determine mutation spectrum
- **Example:** `uniques_only: True`

</details>

<details>
<summary>NanoPlot Integration</summary>

Parameters for NanoPlot quality control visualization.

**nanoplot**
- **Type:** Boolean
- **Default:** `False`
- **Description:** Enable NanoPlot for sequence quality visualization
- **Example:** `nanoplot: True`

**nanoplot_flags**
- **Type:** String
- **Default:** `'--plots dot'`
- **Description:** Command line flags for NanoPlot. `-o` and `-p` are already added
- **Example:** `nanoplot_flags: '--plots kde --format png'`

</details>

## Configuration File Examples

<details>
<summary>Basic Configuration Example</summary>

```yaml
# Minimal configuration for simple analysis
runs: tags.csv
metadata: metadata
sequences_dir: data
fastq_dir: fastq_pass

# Basic plotting
plot_demux: True
plot_mutation-distribution: True
plot_mutation-frequencies: True

# Thread allocation
threads_alignment: 4
threads_demux: 2
```

</details>

<details>
<summary>Advanced Configuration Example</summary>

```yaml
# Advanced configuration with consensus and time series
runs: tags.csv
timepoints: timepoints.csv
metadata: metadata
sequences_dir: /path/to/central/sequences
fastq_dir: fastq_pass, fastq_fail

# UMI consensus
UMI_mismatches: 2
UMI_consensus_minimum: 3
UMI_consensus_maximum: 10

# RCA consensus
peak_finder_settings: '23,3,27,2'
RCA_consensus_minimum: 3
RCA_consensus_maximum: 20

# Demultiplexing
demux_screen_failures: True
demux_threshold: 0.05

# Analysis parameters
mutation_analysis_quality_score_minimum: 10
sequence_length_threshold: 0.1

# Plotting
colormap: 'fire'
plot_genotypes2D: True
genotypes2D_plot_all: True
export_SVG: True

# Dashboard
dashboard_input: experiment1
dashboard_port: 8080

# Performance
threads_alignment: 8
threads_demux: 4
threads_medaka: 4
```

</details>