## v0.8.6
# Major
- Added an example fastq.gz file
- removed support for Guppy basecaller (difficult to maintain and live basecalling is much preferred)
# Minor
- Batch size option for RCA in config
- Providing a dictionary of tag:flags to minimap flags allows for using different minimap settings for different tags
- Automated decision to merge paired end reads or combine batches of reads from a folder (i.e. a nanopore sequencing run folder)