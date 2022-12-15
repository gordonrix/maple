## v0.8.9
# Minor
 - Fixed genotype plotting color/sizes

## v0.8.8
# Major
 - Plotting of demux counts
 - Reorganization of Snakemake modules

## v0.8.7
# Major
- Plotting of RCA consensus stats
- Sequence clustering / dimensionality reduction using PaCMAP
- Scatter plotting of dimension-reduced genotypes

## v0.8.6
# Major
- Added an example fastq.gz file
- removed support for Guppy basecaller (difficult to maintain and live basecalling is much preferred)
# Minor
- Batch size option for RCA in config
- Providing a dictionary of tag:flags to minimap flags allows for using different minimap settings for different tags
- Automated decision to merge paired end reads or combine batches of reads from a folder (i.e. a nanopore sequencing run folder)