## v0.9.2
# Minor
 - Fixed dashboard bug resulting in incorrect muts_of_interest in dashboard
 - Dramatic speedup of muts_of_interest calculation in dashboard

## v0.9.1
# Minor
 - Fixed dashboard bug resulting in incorrect selections when downsampling

## v0.9.0
# Major
 - Added a holoviz-based interactive dashboard
# Minor
 - Added amino acid level sequence embeddings for dimensionality reduction (DR)
 - Fixed how sequence embeddings are done for DR. Mutations are now encoded in binary format.

## v0.8.11
# Minor
 - Simplified sequence import

## v0.8.10
# Minor
 - Fixed color mapping for RCA plot
 - Added timepoint and numerical color mapping for genotypes2D
 - Refactored and improved plotting of distributions of hamming distance / mutations / demux counts / UMI subreads
 - Added hamming distance and mutation distribution plots for timepoints

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