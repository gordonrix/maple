## v0.11.0
# Major
   - Support for MacOS (ARM)
# Minor
   - Added genotypes distribution plot
   - Added support for .csv input for tags and timepoints
   - Bugfixes for empty and nearly empty sequence datasets (demux.py)
   - Default folder that holds metadata changed from 'ref' to 'metadata'
   - Timepoints now defined globally, not for each tag
   - Notices that can be repetitive are now sent to the log file
   - Remove requirement for nanopore option in config
   - Remove references to Singularity

## v0.10.3
# Minor
   - Changed how timepoint files are requested to prevent a bug
   - Fixed a bug in SequenceAnalyzer.downsample()
   - Updated to Python 3.9
   - Added sample/file specificity to SVG export
   - Enrichment now works for timepoint samples that use more than one tag
   - Color by barcode in dashboard
   - Dashboard bug fixes
   - Conda environment: update panel, add jupyter

## v0.10.2
# Major
   - Implemented enrichment scores for noSplit barcodes. A plot comparing replicates will be automatically
      generated, and when NT or AA analysis are performed,
      scores are appended to the genotypes CSV and become viewable within the dashboard
# Minor
   - Updated README
   - Enabled usage without sequence import (sequences must be named the same as the tags and in the sequences folder)
   - Custom colormaps using config input
   - Dashboard:
      - can now color and filter by any numerical columns
      - histogram now works for any numerical column
      - points plot axes can now be any numerical column
      - enrichment scores are now integrated into the dashboard
         if they are calculated for the sample being used for the dashboard

## v0.10.1
# Minor
 - Dashboard:
    - improved consensus export
    - terminal feedback on export

## v0.10.0
# Major
 - introduced SequenceAnalyzer class for interconversion of multiple sequences
    among different encodings and fast, vectorized analysis operations
 - New dashboard features:
    - implemented SequenceAnalyzer backend
    - datashaded points plot (now default)
    - data export (SVG plots, CSVs, sequences as fasta, consensus sequences from selection)
# Minor
 - New required package spatialpandas, for making complex selections with datashader
 - Mutation analysis genotypes output now includes AA insertions and deletions

## v0.9.3
# Major
 - Reworked mutation frequency plotting, including most/least common mutations
 - Added violin mutation distribution plots
 - Enable SVG export for most plots 
# Minor
 - Genotypes 2D hexbins plot

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