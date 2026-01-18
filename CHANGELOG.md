## v1.0.1
# Minor
   - Fixed genotype enrichment pipeline
   - Fixed demux enrichment merging with genotypes
   - Added GitHub Actions CI workflow for Linux and macOS testing
   - Documentation updates

## v1.0.0
# Major
   - **Multi-reference support**: Support for analyzing multiple reference sequences simultaneously
      - Multi-reference demultiplexing, consensus generation, and mutation analysis
      - Reference-specific visualization with stacked bar charts and gridded distribution plots
      - CSV-based reference configuration system
   - **Enhanced enrichment pipeline**: Major improvements to enrichment analysis
      - Genotype enrichment analysis with filtering and validation
      - Multi-reference enrichment support
      - Full integration with dashboard for interactive exploration
      - Improved checkpoint handling and dimension reduction
   - **Dashboard overhaul**: Complete rebuild with more robust Panel reactive architecture
      - Multi-dataset support for comparing multiple experiments
      - PDB/CIF structure viewer integration with py3Dmol
      - Improved performance and user experience
      - Enhanced mutation visualization and filtering
   - **UMI workflow refactor**: Consensus pipeline modernized to use BAM format
      - New UMI_process.py for efficient multi-reference consensus generation
      - Improved memory usage and processing speed
      - Better handling of demuxed BAM files
   - **Streamlined installation**: New install folder with lock files for reproducible environments
   - **Configuration improvements**: Enhanced validation and better organization
      - Timepoint column validation for enrichment analysis
      - Improved CSV format validation with backwards compatibility
      - Better config file management (config_final.yaml)
   - **Modularized custom scripts**: The following scripts have been reworked to enable standalone use:
      - mutation_analysis.py
      - mutation_statistics.py
      - demux.py
      - enrichment.py
      - UMI_process.py
      - UMI_groups_log.py
      - generate_barcode_ref.py
      - plot_demux.py
      - plot_distribution.py
      - structure_heatmap.py

# Minor
   - Memory optimization: SequenceAnalyzer no longer stores full one-hot encoding on initialization
   - Performance improvements:
      - RCA consensus now uses an optimized fork of C3POa.py
      - UMI consensus avoids aligning twice
      - Optimized non-exact barcode matching performance
      - Better handling of large datasets in dashboard
   - Better error handling:
      - Empty output handling for mutation analysis
      - Fix for DtypeWarning in pandas CSV reads
      - Improved validation throughout pipeline
   - Documentation improvements:
      - Enhanced timepoints CSV documentation
      - Updated configuration documentation
      - Better example metadata and working directory
   - Barcode generation and validation:
      - Validation for generated barcode consistency across tags
      - Improved barcode tracking system
   - Utility improvements:
      - New compile_original_reads utility
      - Enhanced enrichment utility functions
      - Removed deprecated utility scripts
   - Various bugfixes

## v0.11.0
# Major
   - Support for MacOS (ARM)
   - Switched to UMIcollapse for UMI processing (large speedup for big UMI datasets)
   - barcodeGroups converted to partition_barcode_groups, added label_barcode groups
      - label_barcode_groups use barcode groups to label demux output sequences
      - Can now assign multiple barcode groups to the same partition/label barcode group
   - demux.py refactor
# Minor
   - Added genotypes distribution plot
   - Added support for .csv input for tags and timepoints
   - Bugfixes for empty and nearly empty sequence datasets (demux.py)
   - Default folder that holds metadata changed from 'ref' to 'metadata'
   - Timepoints now defined globally, not for each tag
   - Notices that can be repetitive are now sent to the log file
   - Remove requirement for nanopore option in config
   - Remove references to Singularity
   - Separate barcode tracking from genotype tracking (currently, breaks merge_enrichment.py)
   - Fixed bug that would cause the 'representative alignment' .txt output from mutation_analysis.py to be incorrect
   - Fixed aesthetic issues with SVG export
   - Updated package versions to install
   - Updated to latest Snakemake version, adjusted targets and mutation_analysis rules to accomodate this

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