"""Creation of consensus sequences from repetitive reads.
Differs from smolecule.py from medaka in that it uses BAM files as input,
taking the reference sequence as the initial consensus instead of generating one de novo,
and using the alignments directly instead of generating new alignments.
Compatible with medaka v2.1.1.
"""
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor
import functools
import h5py
import logging
import os
import re
import time
from timeit import default_timer as now
import warnings

import mappy
import numpy as np
import parasail
import pysam
import spoa

import medaka.align
import medaka.common
import medaka.medaka

Subread = namedtuple('Subread', 'name seq')
Alignment = namedtuple('Alignment', 'rname qname flag rstart seq cigar')


def write_bam(fname, alignments, header, bam=True):
    """Write a `.bam` file for a set of alignments.

    :param fname: output filename.
    :param alignments: a list of `Alignment` tuples.
    :param header: bam header
    :param bam: write bam, else sam

    """
    mode = 'wb' if bam else 'w'
    if isinstance(header, dict):
        header = pysam.AlignmentHeader.from_dict(header)
    with pysam.AlignmentFile(fname, mode, header=header) as fh:
        for subreads in alignments:
            for aln in sorted(subreads, key=lambda x: x.rstart):
                a = medaka.align.initialise_alignment(
                    aln.qname, header.get_tid(aln.rname),
                    aln.rstart, aln.seq,
                    aln.cigar, aln.flag)
                fh.write(a)
    if mode == 'wb':
        pysam.index(fname)


def count_processed_samples(hdf_file):
    """Count number of samples/chunks processed in medaka HDF5 file.

    :param hdf_file: path to medaka consensus HDF5 file.
    :returns: number of samples written, or 0 if file doesn't exist or can't be read.
    """
    try:
        with h5py.File(hdf_file, 'r') as f:
            # Count top-level groups (each consensus is stored as a group)
            # Try multiple approaches to count samples
            if 'sample_registry' in f:
                # sample_registry exists - count entries
                registry = f['sample_registry']
                return len(registry) if hasattr(registry, '__len__') else 0
            else:
                # Count top-level groups
                return len([k for k in f.keys() if not k.startswith('_')])
    except Exception:
        return 0


def read_bam_workflow(bam_path, references):
    """Read BAM file and extract alignment information, bypassing mappy alignment.

    This function reads a BAM file produced by UMI_process with alignments already
    present, groups reads by UMI, and prepares them for medaka consensus without
    re-aligning.

    Read names in BAM must follow format: UMI-{unique_id}|ref-{ref_name}|{read_id}

    :param bam_path: Path to input BAM file with aligned reads
    :param references: Dict of reference sequences {ref_name: sequence}
    :returns: header, consensuses, alignments (same format as poa_workflow)
    """
    logger = medaka.common.get_named_logger('BAMReader')

    # Parse BAM file and group alignments by (UMI, reference)
    logger.info(f"Reading BAM file: {bam_path}")
    bam_in = pysam.AlignmentFile(bam_path, 'rb')

    # Group alignments: {(ref_name, umi_id): [AlignedSegment, ...]}
    umi_alignments = {}

    for read in bam_in.fetch():
        if read.is_unmapped:
            continue

        # Parse read name: UMI-{id}|ref-{ref_name}|{read_id}
        try:
            parts = read.query_name.split('|')
            umi_part = parts[0]  # UMI-{id}
            ref_part = parts[1]  # ref-{ref_name}

            umi_id = umi_part.split('-', 1)[1]  # Extract {id}
            ref_name = ref_part.split('-', 1)[1]  # Extract {ref_name}

            key = (ref_name, umi_id)
            if key not in umi_alignments:
                umi_alignments[key] = []
            umi_alignments[key].append(read)

        except (IndexError, ValueError) as e:
            logger.warning(f"Could not parse read name: {read.query_name}")
            continue

    bam_in.close()

    logger.info(f"Found {len(umi_alignments)} UMI groups from BAM")

    # Build header and prepare outputs in same format as poa_workflow
    header = {'HD': {'VN': 1.0}, 'SQ': []}
    consensuses = []
    alignments = []

    for (ref_name, umi_id), bam_reads in umi_alignments.items():
        # Get reference sequence
        if ref_name not in references:
            logger.warning(f"Reference '{ref_name}' not found, skipping UMI group {umi_id}")
            continue

        consensus_seq = references[ref_name]
        rname = f"UMI-{umi_id}|ref-{ref_name}"

        # Add to header
        header['SQ'].append({
            'LN': len(consensus_seq),
            'SN': rname
        })

        # Add consensus (using reference as consensus since no POA)
        consensuses.append([rname, consensus_seq])

        # Convert BAM reads to Alignment namedtuples
        aligns_for_group = []
        for bam_read in bam_reads:
            # Convert pysam AlignedSegment to medaka Alignment format
            aln = Alignment(
                rname=rname,
                qname=bam_read.query_name,
                flag=bam_read.flag,
                rstart=bam_read.reference_start,
                seq=bam_read.query_sequence,
                cigar=bam_read.cigarstring
            )
            aligns_for_group.append(aln)

        alignments.append(aligns_for_group)

    logger.info(
        "Extracted {} consensus with {} alignments from BAM.".format(
            len(consensuses), len(alignments)))

    return header, consensuses, alignments


class MyArgs:
    """Wrap the given args with defaults for prediction function."""

    # these need to be class attrs
    args = None
    defaults = None

    def __init__(self, args, defaults):
        """Initialize the class."""
        self.args = args
        self.defaults = defaults

    def __getattr__(self, attr):
        """Get the standard argument, or default if not present."""
        try:
            return getattr(self.args, attr)
        except AttributeError:
            return getattr(self.defaults, attr)


def run_medaka_predict(args, output_dir):
    """Run medaka.prediction.predict with stderr redirected to log file.

    This must be a module-level function (not nested) to be pickleable for multiprocessing.
    Redirects stderr to suppress warnings like "Pileup counts do not span requested region".
    """
    import sys
    medaka_log = os.path.join(output_dir, 'medaka.log')
    with open(medaka_log, 'w') as log_fh:
        old_stderr = sys.stderr
        sys.stderr = log_fh
        try:
            result = medaka.prediction.predict(args)
        finally:
            sys.stderr = old_stderr
    return result


def main(args):
    """Entry point for repeat read consensus creation."""
    parser = medaka.medaka.medaka_parser()
    defaults = parser.parse_args([
        "inference", medaka.medaka.CheckBam.fake_sentinel,
        "fake_out"])

    args = MyArgs(args, defaults)

    # Setup logging: INFO to console with batch identifier, all levels to file
    batch_name = os.path.basename(args.output)

    # Create output directory if it doesn't exist
    medaka.common.mkdir_p(args.output, info='Results will be overwritten.')

    # File handler - captures everything
    log_file = os.path.join(args.output, 'maple_smolecule.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))

    # Console handler - only INFO and above, with batch name prefix
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter(f'[{batch_name}] %(message)s'))

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    root_logger.handlers = []  # Clear any existing handlers
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)

    # Get maple_smolecule logger
    logger = medaka.common.get_named_logger('Smolecule')
    logger.setLevel(logging.INFO)

    if args.chunk_ovlp >= args.chunk_len:
        raise ValueError(
            'chunk_ovlp {} must be smaller than chunk_len {}'.format(
                args.chunk_ovlp, args.chunk_len))

    # Read reference(s) from alignment file into dict
    references = {}
    with pysam.FastxFile(args.reference) as fh:
        for entry in fh:
            references[entry.name] = entry.sequence.upper()

    logger.debug("Loaded {} reference(s).".format(len(references)))

    # Use BAM workflow - read alignments directly
    input_file = args.fasta[0] if args.fasta else None
    logger.info("Using BAM input with pre-existing alignments.")

    t0 = now()
    header, consensuses, alignments = read_bam_workflow(
        input_file,
        references
    )
    t1 = now()

    logger.info(
        "Writing medaka input bam for {} reads.".format(len(alignments)))
    bam_file = os.path.join(args.output, 'subreads_to_spoa.bam')
    write_bam(bam_file, alignments, header)

    spoa_file = os.path.join(args.output, 'poa.fasta')
    with open(spoa_file, 'w') as fh:
        for rname, cons in consensuses:
            fh.write('>{}\n{}\n'.format(rname, cons))

    total_samples = len(consensuses)
    logger.info(f"Running medaka consensus on {total_samples} reads.")
    t2 = now()
    args.bam = bam_file
    out_dir = args.output
    args.output = os.path.join(out_dir, 'consensus.hdf')

    # Run medaka prediction in a subprocess so GPU resources are cleaned up
    with ProcessPoolExecutor() as executor:
        fut = executor.submit(run_medaka_predict, args, out_dir)
        _ = fut.result()
    t3 = now()

    logger.info("Running medaka stitch.")
    args.draft = spoa_file
    args.inputs = [args.output]
    out_ext = 'fasta'
    if args.qualities:
        out_ext = 'fastq'
    args.output = os.path.join(out_dir, 'consensus.{}'.format(out_ext))
    args.regions = None  # medaka consensus changes args.regions
    args.fillgaps = False
    # NOTE: medaka.stitch.stitch() appends '_0' suffix to all sequence IDs in the output
    medaka.stitch.stitch(args)
    logger.info(
        "Single-molecule consensus sequences written to {}.".format(
            args.output))
    logger.info("Medaka run time: {:.0f}s".format(t3 - t2))
