import argparse
import logging
import os

import pysam

import medaka.common
import medaka.datastore
import medaka.features
import medaka.labels
import medaka.models
import medaka.options
import medaka.prediction
import medaka.rle
import medaka.maple_smolecule
import medaka.smolecule
import medaka.stitch
import medaka.training
import medaka.variant
import medaka.vcf
import medaka.wrappers


class ResolveModel(argparse.Action):
    """Resolve model filename or ID into filename"""
    def __init__(
            self, option_strings, dest, default=None, required=False,
            help='Model file.'):
        super().__init__(
            option_strings, dest, nargs=1, default=default, required=required,
            help='{} {{{}}}'.format(help, ', '.join(medaka.options.allowed_models)))

    def __call__(self, parser, namespace, values, option_string=None):
        val = values[0]
        try:
            model_fp = medaka.models.resolve_model(val)
        except Exception as e:
            msg = "Error validating model from '--{}' argument: {}."
            raise RuntimeError(msg.format(self.dest, str(e)))
        #TODO: verify the file is a model?
        setattr(namespace, self.dest, model_fp)


class CheckBlockSize(argparse.Action):
    """Check that block_size < 94"""

    def __call__(self, parser, namespace, values, option_string=None):
        if values > 94:
            parser.error(
                'Maximum block_size is 94, to avoid going over the ASCII '
                'limit of 127 (scores start in ASCII character 33).')
        setattr(namespace, self.dest, values)


class CheckBam(argparse.Action):
    """Check a bam is a bam."""
    fake_sentinel = "CheckBAMFake.bam"  # see below

    def __call__(self, parser, namespace, values, option_string=None):
        # allow us to skip the file exist check. This is used to allow
        # smolecule to run the consensus argument parser to inherit arguments
        # required to run the prediction program without itself having all
        # the arguments on its CLI
        if values == self.fake_sentinel:
            setattr(namespace, self.dest, values)
            return

        if not os.path.exists(values):
            raise RuntimeError(
                "Filepath for '--{}' argument does not exist ({})".format(
                    self.dest, values)
            )
        try:
            with pysam.AlignmentFile(values) as bam:
                # just to check its a bam
                _ = bam.references
        except Exception:
            raise IOError('The bam {} could not be read.'.format(values))
        setattr(namespace, self.dest, values)

    @staticmethod
    def check_read_groups(fname, rg=None):
        """Check bam read groups are consistent with the specified read group.

        Raises a RuntimeError if:
            * no read group was specified but the bam has read groups
            * user specified a read group but the bam has no read groups
            * user specified a read group but it it not amongst bam read groups

        :param fname: bam file name.
        :param rg: check if this read group is present.

        :raises: RuntimeError
        """
        with pysam.AlignmentFile(fname) as bam:
            # As of 13/12/19 pypi still has no wheel for pysam v0.15.3 so we
            # pinned to v0.15.2. However bioconda's v0.15.2 package
            # conflicts with the libdeflate they have so we are forced
            # to use newer versions. Newer versions however have a
            # different API, sigh...
            try:
                header_dict = bam.header.as_dict()
            except AttributeError:
                header_dict = bam.header

            if 'RG' in header_dict:
                read_groups = set([r['ID'] for r in header_dict['RG']])
            else:
                read_groups = set()
            if rg is not None and len(read_groups) == 0:
                # User asked for a read group but none were found
                raise RuntimeError('No RG tags found in the bam {}'.format(fname))
            elif rg is None and len(read_groups) > 1:
                # User did not ask for a read group but groups are present.
                raise RuntimeError(
                    'The bam {} contains more than one read group. '
                    'Please specify `--RG` to select which read group'
                    'to process from {}'.format(fname, read_groups))
            elif rg is not None and len(read_groups) > 0:
                # User asked for a read group and groups are present.
                if rg not in read_groups:
                    msg = 'RG {} is not in the bam {}. Try one of {}'
                    raise RuntimeError(msg.format(rg, fname, read_groups))


class RegionParser(argparse.Action):
    """Parse regions, checking if --regions option is a bed file."""

    def __call__(self, parser, namespace, values, option_string=None):
        if len(values) == 1 and os.path.exists(values[0]):
            # parse bed file
            regions = []
            for chrom, start, stop in medaka.common.yield_from_bed(values[0]):
                regions.append('{}:{}-{}'.format(chrom, start, stop))
        else:
            regions = values
        regions = [medaka.common.Region.from_string(r) for r in regions]
        setattr(namespace, self.dest, regions)


class RegionRefNameParser(RegionParser):
    """Parse regions, retaining only region ref_names."""

    def __call__(self, parser, namespace, values, option_string=None):
        super().__call__(parser, namespace, values, option_string=None)
        regions = getattr(namespace, self.dest)
        if any((r.start is not None or r.end is not None for r in regions)):
            print('WARNING: This program can only process entire ' +
                  'contigs, ignoring region start and end coordinates.')
        # keep regions in order but avoid duplicating ref_names.
        ref_names = []
        for r in regions:
            if r.ref_name not in ref_names:
                ref_names.append(r.ref_name)
        regions = [medaka.common.Region(r, None, None) for r in ref_names]
        setattr(namespace, self.dest, regions)


def _log_level():
    """Parser to set logging level and acquire software version/commit"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

    #parser.add_argument('--version', action='version', version=get_version())

    modify_log_level = parser.add_mutually_exclusive_group()
    modify_log_level.add_argument('--debug', action='store_const',
        dest='log_level', const=logging.DEBUG, default=logging.INFO,
        help='Verbose logging of debug information.')
    modify_log_level.add_argument('--quiet', action='store_const',
        dest='log_level', const=logging.WARNING, default=logging.INFO,
        help='Minimal logging; warnings only).')

    return parser


def _model_arg():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--model', action=ResolveModel,
            default=medaka.options.default_models['consensus'],
            help='Model to use.')
    return parser


def _min_depth_arg():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--min_depth', type=int, default=0,
        help="Sites with depth lower than min_depth will not be polished.")
    return parser


def _rg_arg():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    rg_group = parser.add_argument_group('read group', 'Filtering alignments the read group (RG) tag, expected to be string value.')
    rg_group.add_argument('--RG', metavar='READGROUP', type=str, help='Read group to select.')
    return parser


def _align_chunking():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    group = parser.add_argument_group('CONSENSUS2VCF', 'consensus to '
                                      + ' reference alignment options.')
    group.add_argument('--chunk_size', type=int, default=100000,
                       help='Size of consensus chunks')
    group.add_argument('--pad', type=int, default=10000,
                       help='Reference chunks are chunk_size + pad.')
    group.add_argument('--mode', default='HWT', choices=['NW', 'HW', 'HWT'],
                       help='Edlib alignment mode. ' +
                       'NW: global in consensus and ref. ' +
                       'HW: global in consensus, local in ref. ' +
                       'HWT: same as HW, but alignments trimmed to ' +
                       'start and end on a match.'
                       )
    return parser


def _regions_or_bed_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--regions', default=None, action=RegionParser, nargs='+',
                        help='Genomic regions to analyse, or a bed file.')
    return parser


def _region_ref_names():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--regions', default=None, action=RegionRefNameParser,
                        nargs='+',
                        help='Genomic ref_names to process, or a bed file.')
    return parser


def _chunking_feature_args(batch_size=100, chunk_len=10000, chunk_ovlp=1000):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)
    parser.add_argument('--batch_size', type=int, default=batch_size, help='Inference batch size.')
    parser.add_argument('--chunk_len', type=int, default=chunk_len, help='Chunk length of samples.')
    parser.add_argument('--chunk_ovlp', type=int, default=chunk_ovlp, help='Overlap of chunks.')
    return parser


def _validate_common_args(args):
    """Do some common argument validation."""
    logger = medaka.common.get_named_logger('ValidArgs')
    if hasattr(args, 'bam') and args.bam is not None:
        RG = args.RG if hasattr(args, 'RG') else None
        CheckBam.check_read_groups(args.bam, RG)
        if RG is not None:
            msg = "Reads will be filtered to only those with RG tag: {}"
            logger.info(msg.format(RG))
    # if model is default, resolve to file, save mess in help text
    if hasattr(args, 'model') and args.model is not None:
        args.model = medaka.models.resolve_model(args.model)


def print_model_path(args):
    print(os.path.abspath(args.model))


def is_rle_model(args):
    print(is_rle_encoder(args.model))


def is_rle_encoder(model_name):
    """ Return encoder used by model"""
    rle_encoders = [medaka.features.HardRLEFeatureEncoder]
    modelstore = medaka.models.open_model(model_name)
    encoder = modelstore.get_meta('feature_encoder')
    is_rle = issubclass(type(encoder), medaka.features.HardRLEFeatureEncoder)
    return is_rle


def get_alignment_params(model):
    if is_rle_encoder(model):
        align_params = medaka.options.alignment_params['rle']
    else:
        align_params = medaka.options.alignment_params['non-rle']
    return align_params


def print_alignment_params(args):
    print(get_alignment_params(args.model))


def print_all_models(args):
    print('Available:', ', '.join(medaka.options.allowed_models))
    for key in ('consensus', 'variant'):
        # medaka_variant relies on this order
        print('Default {}: '.format(key), medaka.options.default_models[key])


def fastrle(args):
    import libmedaka
    libmedaka.lib.fastrle(args.input.encode(), args.block_size)


def download_models(args):
    logger = medaka.common.get_named_logger('ResolveModels')
    for model in args.models:
        ns = argparse.Namespace()
        fp = ResolveModel('--model', 'model')(None, ns, [model], 'model')
        logger.info('Model {} is stored at: {}'.format(model, ns.model))


class StoreDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        """Converts 'key1=value1 key2=value1a,value2a' to a dictionary.

        List values are supported, as are simple type conversions.
        """
        str_to_type = {
            'None': None,
            'True': True, 'False': False,
            'true': True, 'false': False,
            'TRUE': True, 'FALSE': False}

        def _str_to_numeric(x):
            if not isinstance(x, str):
                return x
            try:
                return int(x)
            except:
                try:
                    return float(x)
                except:
                    return x

        my_dict = {}
        for kv in values:
            key, value = kv.split("=")
            list_like = ',' in value
            value = str_to_type.get(value, value)
            if value is not None:
                if list_like:
                    value = [_str_to_numeric(str_to_type.get(x,x))
                        for x in value.split(',')]
                else:
                    value = _str_to_numeric(value)
            my_dict[key] = value
        setattr(namespace, self.dest, my_dict)


def medaka_parser():
    """Create the medaka command-line interface."""
    from medaka import __version__
    parser = argparse.ArgumentParser('medaka',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid commands', help='additional help', dest='command')
    subparsers.required = True

    parser.add_argument('--version', action='version',
        version='%(prog)s {}'.format(__version__))

    # Compress bam file_ext
    rparser = subparsers.add_parser('compress_bam',
        help='Compress an alignment into RLE form. ',
        parents=[_log_level(), _regions_or_bed_args()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rparser.set_defaults(func=medaka.rle.compress_bam)
    rparser.add_argument('bam_input', help='Bam file to compress.')
    rparser.add_argument('bam_output', help='Output bam file.')
    rparser.add_argument('ref_fname',
                         help='Reference fasta file used for `bam_input`.')
    rparser.add_argument('--threads', type=int, default=1,
                         help='Number of threads for parallel execution.')

    rparser.add_argument(
        '--use_fast5_info', metavar='<fast5_dir> <index>', default=None,
        nargs=2, help=(
            'Root directory containing the fast5 files and .tsv file with '
            'columns filename and read_id.'))

    # Creation of feature files
    fparser = subparsers.add_parser('features',
        help='Create features for inference.',
        parents=[_log_level(), _chunking_feature_args(), _regions_or_bed_args()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    fparser.set_defaults(func=medaka.features.create_samples)
    fparser.add_argument('bam', help='Input alignments.', action=CheckBam)
    fparser.add_argument('output', help='Output features file.')
    fparser.add_argument('--truth', help='Bam of truth aligned to ref to create features for training.')
    fparser.add_argument('--truth_haplotag', help='Two-letter tag defining haplotype of alignments for polyploidy labels.')
    fparser.add_argument('--threads', type=int, default=1, help='Number of threads for parallel execution.')
    # TODO: enable other label schemes.
    fparser.add_argument('--label_scheme', default='HaploidLabelScheme', help='Labelling scheme.',
                         choices=sorted(medaka.labels.label_schemes))
    fparser.add_argument('--label_scheme_args', action=StoreDict, nargs='+',
        default=dict(), metavar="KEY1=VAL1 KEY2=VAL2a,VAL2b...",
        help="Label scheme key-word arguments.")
    fparser.add_argument('--feature_encoder', default='CountsFeatureEncoder',
        help='FeatureEncoder used to create the features.',
        choices=sorted(medaka.features.feature_encoders))
    fparser.add_argument('--feature_encoder_args', action=StoreDict, nargs='+',
        default=dict(), metavar="KEY1=VAL1 KEY2=VAL2a,VAL2b...",
        help="Feature encoder key-word arguments.")

    # Training program
    tparser = subparsers.add_parser('train',
        help='Train a model from features.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    tparser.set_defaults(func=medaka.training.train)
    tparser.add_argument('features', nargs='+', help='Paths to training data.')
    tparser.add_argument('--train_name', type=str, default='keras_train', help='Name for training run.')
    tparser.add_argument('--model', action=ResolveModel, help='Model definition and initial weights .hdf, or .yml with kwargs to build model.')
    tparser.add_argument('--epochs', type=int, default=5000, help='Maximum number of trainig epochs.')
    tparser.add_argument('--batch_size', type=int, default=100, help='Training batch size.')
    tparser.add_argument('--max_samples', type=int, default=float("inf"), help='Only train on max_samples.')
    tparser.add_argument('--mini_epochs', type=int, default=1, help='Reduce fraction of data per epoch by this factor')
    tparser.add_argument('--seed', type=int, help='Seed for random batch shuffling.')
    tparser.add_argument('--threads_io', type=int, default=1, help='Number of threads for parallel IO.')
    tparser.add_argument('--device', type=int, default=0, help='GPU device to use.')
    tparser.add_argument('--optimizer', type=str, default='rmsprop', choices=['nadam','rmsprop'], help='Optimizer to use.')
    tparser.add_argument('--optim_args', action=StoreDict, default=None, nargs='+',
        metavar="KEY1=VAL1,KEY2=VAL2...", help="Optimizer key-word arguments.")

    vgrp = tparser.add_mutually_exclusive_group()
    vgrp.add_argument('--validation_split', type=float, default=0.2, help='Fraction of data to validate on.')
    vgrp.add_argument('--validation_features', nargs='+', default=None, help='Paths to validation data')

    # Consensus from bam input
    # NB: any args added here should be set to default values in smolecule:main()
    #     to avoid attribute errors in that program.
    cparser = subparsers.add_parser('consensus',
        help='Run inference from a trained model and alignments.',
        parents=[_log_level(), _chunking_feature_args(), _regions_or_bed_args(), _model_arg(), _rg_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    cparser.set_defaults(func=medaka.prediction.predict)
    cparser.add_argument('bam', help='Input alignments.', action=CheckBam)
    cparser.add_argument('output', help='Output file.')
    cparser.add_argument('--threads', type=int, default=1, help='Number of threads used by inference.')
    cparser.add_argument('--bam_workers', type=int, default=2, help='Number of workers used to prepare data from bam.')
    cparser.add_argument('--bam_chunk', type=int, default=int(1e6), help='Size of reference chunks each worker parses from bam. (can be used to control memory use).')
    cparser.add_argument('--check_output', action='store_true', default=False,
            help='Verify integrity of output file after inference.')
    cparser.add_argument('--save_features', action='store_true', default=False,
            help='Save features with consensus probabilities.')
    tag_group = cparser.add_argument_group('filter tag', 'Filtering alignments by an integer valued tag.')
    tag_group.add_argument('--tag_name', type=str, help='Two-letter tag name.')
    tag_group.add_argument('--tag_value', type=int, help='Value of tag.')
    tag_group.add_argument('--tag_keep_missing', action='store_true', help='Keep alignments when tag is missing.')

    # Consensus from single-molecules with subreads, BAM file as input
    msparser = subparsers.add_parser('maple_smolecule',
        help='Create consensus sequences from pre-aligned single-molecule reads.',
        parents=[_log_level(), _chunking_feature_args(batch_size=100, chunk_len=1000, chunk_ovlp=500), _model_arg(), _min_depth_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    msparser.set_defaults(func=medaka.maple_smolecule.main)
    msparser.add_argument('output', help='Output directory.')
    msparser.add_argument('reference', help='fasta file for reference sequence that will be used for alignment.')
    msparser.add_argument('fasta', nargs='+', help='pre-aligned single-molecule reads, sorted and grouped by UMI ID (UG tag in BAM file).')
    msparser.add_argument('--method', choices=['spoa'], default='spoa', help='Pre-medaka consensus generation method.')
    msparser.add_argument('--depth', type=int, default=3, help='Minimum subread count.')
    msparser.add_argument('--length', type=int, default=400, help='Minimum median subread length.')
    msparser.add_argument('--threads', type=int, default=1, help='Number of threads used by inference.')
    msparser.add_argument('--check_output', action='store_true', default=False,
            help='Verify integrity of output file after inference.')
    msparser.add_argument('--save_features', action='store_true', default=False,
            help='Save features with consensus probabilities.')
    msparser.add_argument('--qualities', action='store_true', default=False,
            help='Output consensus with per-base quality scores (fastq).')

    # Consensus from single-molecules with subreads
    smparser = subparsers.add_parser('smolecule',
        help='Create consensus sequences from single-molecule reads.',
        parents=[_log_level(), _chunking_feature_args(batch_size=100, chunk_len=1000, chunk_ovlp=500), _model_arg(), _min_depth_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    smparser.set_defaults(func=medaka.smolecule.main)
    smparser.add_argument('output', help='Output directory.')
    smparser.add_argument('fasta', nargs='+', help='Single-molecule reads, one file per read.')
    smparser.add_argument('--method', choices=['spoa'], default='spoa', help='Pre-medaka consensus generation method.')
    smparser.add_argument('--depth', type=int, default=3, help='Minimum subread count.')
    smparser.add_argument('--length', type=int, default=400, help='Minimum median subread length.')
    smparser.add_argument('--threads', type=int, default=1, help='Number of threads used by inference.')
    smparser.add_argument('--check_output', action='store_true', default=False,
            help='Verify integrity of output file after inference.')
    smparser.add_argument('--save_features', action='store_true', default=False,
            help='Save features with consensus probabilities.')
    smparser.add_argument('--qualities', action='store_true', default=False,
            help='Output consensus with per-base quality scores (fastq).')

    # Consensus from features input
    cfparser = subparsers.add_parser('consensus_from_features',
        help='Run inference from a trained model on existing features.',
        parents=[_log_level(), _model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # TODO: why doesn't this have a .set_defaults?
    cfparser.add_argument('features', nargs='+', help='Pregenerated features (from medaka features).')

    # Compression of fasta/q to quality-RLE fastq
    rleparser = subparsers.add_parser('fastrle',
        help='Create run-length encoded fastq (lengths in quality track).',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rleparser.set_defaults(func=fastrle)
    rleparser.add_argument('input', help='Input fasta/q. may be gzip compressed.')
    rleparser.add_argument('--block_size', action=CheckBlockSize, default=94, type=int,
        help='Block size for hompolymer splitting, e.g. with a value of blocksize=3, AAAA -> A3 A1.')

    # Post-processing of consensus outputs
    sparser = subparsers.add_parser('stitch',
        help='Stitch together output from medaka consensus into final output.',
        parents=[_log_level(), _region_ref_names(), _min_depth_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sparser.set_defaults(func=medaka.stitch.stitch)
    sparser.add_argument('inputs',
        help='Consensus .hdf files.', nargs='+')
    sparser.add_argument('draft',
        help='Draft .fasta. Consensus gaps will be filled with unpolished '
        'draft sequence to avoid contig fragmentation.')
    sparser.add_argument('output',
        help='Output .fasta.', default='consensus.fasta')
    sparser.add_argument('--threads',
        help='Number of worker processes to use.', default=1, type=int)
    sparser.add_argument('--no-fillgaps',
        help="Don't fill gaps in consensus sequence with draft sequence.",
        default=True, action='store_false', dest='fillgaps')
    sparser.add_argument(
        '--fill_char',
        help="Use a designated character to fill gaps.",
        default=None, type=str)
    sparser.add_argument('--qualities',
        help="Output with per-base quality scores (fastq).",
        action='store_true')

    var_parser = subparsers.add_parser('variant',
        help='Decode probabilities to VCF.',
        parents=[_log_level(), _region_ref_names()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    var_parser.set_defaults(func=medaka.variant.variants_from_hdf)
    var_parser.add_argument('ref_fasta', help='Reference sequence .fasta file.')
    var_parser.add_argument('inputs', nargs='+', help='Consensus .hdf files.')
    var_parser.add_argument('output', help='Output .vcf.', default='medaka.vcf')
    var_parser.add_argument('--verbose', action='store_true',
                            help='Populate VCF info fields.')
    var_parser.add_argument('--ambig_ref', action='store_true',
                         help='Decode variants at ambiguous reference positions.')
    var_parser.add_argument('--gvcf', action='store_true',
                         help='Output VCF records for reference loci predicted to be non-variant.')

    # TODO do we still need this?
    snp_parser = subparsers.add_parser('snp',
        help='Decode probabilities to SNPs.',
        parents=[_log_level(), _region_ref_names()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    snp_parser.set_defaults(func=medaka.variant.snps_from_hdf)
    snp_parser.add_argument('ref_fasta', help='Reference sequence .fasta file.')
    snp_parser.add_argument('inputs', nargs='+', help='Consensus .hdf files.')
    snp_parser.add_argument('output', help='Output .vcf.', default='medaka.vcf')
    snp_parser.add_argument('--threshold', default=0.04, type=float,
                            help='Threshold for considering secondary calls. A value of 1 will result in haploid decoding.')
    snp_parser.add_argument('--ref_vcf', default=None, help='Reference vcf.')
    snp_parser.add_argument('--verbose', action='store_true',
                         help='Populate VCF info fields.')

    # Tools
    toolparser = subparsers.add_parser('tools',
        help='tools subcommand.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    toolsubparsers = toolparser.add_subparsers(title='tools', description='valid tool commands', help='additional help', dest='tool_command')

    # merge two haploid VCFs into a diploid VCF.
    h2dparser = toolsubparsers.add_parser('haploid2diploid',
        help='Merge two haploid VCFs into a diploid VCF.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    h2dparser.set_defaults(func=medaka.vcf.haploid2diploid)
    h2dparser.add_argument('vcf1', help='Input .vcf file.')
    h2dparser.add_argument('vcf2', help='Input .vcf file.')
    h2dparser.add_argument('ref_fasta', help='Reference .fasta file.')
    h2dparser.add_argument('vcfout', help='Output .vcf.')
    h2dparser.add_argument('--discard_phase', default=False,
                           action='store_true', help='output unphased diploid vcf')
    h2dparser.add_argument('--adjacent', action='store_true',
                         help=('Merge adjacent variants (i.e. variants with non-overlapping genomic ranges) instead' +
                               ' of just overlapping ones. If set to True, all runs of adjacent variants will be' +
                               ' merged including those which appear in just one of the input VCFs.'))
    h2dparser.add_argument('--split_mnp', action='store_true',
                         help='Split each mnp into multiple snps.')

    # split a diploid VCF into two haploid VCFs.
    d2hparser = toolsubparsers.add_parser('diploid2haploid',
        help='Split a diploid VCF into two haploid VCFs.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    d2hparser.set_defaults(func=medaka.vcf.diploid2haploid)
    d2hparser.add_argument('vcf', help='Input .vcf file.')
    d2hparser.add_argument('--notrim', action='store_true',
                         help='Do not trim variant to minimal ref and alt and update pos.')

    # classify variants in a vcf writing new vcfs containing subs, indels etc.
    clparser = toolsubparsers.add_parser('classify_variants',
        help='Classify variants and write vcf for each type and subtype.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    clparser.set_defaults(func=medaka.vcf.classify_variants)
    clparser.add_argument('vcf', help='Input .vcf file.')
    clparser.add_argument('--replace_info', action='store_true',
                         help='Replace info tag (useful for visual inspection of types).')

    # annotate vcf with read depth and supporting reads info
    annparser = toolsubparsers.add_parser('annotate',
        help='Annotate vcf with read depth and supporting reads info fields.',
        parents=[_log_level(), _rg_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    annparser.set_defaults(func=medaka.vcf.annotate_vcf_n_reads)
    annparser.add_argument('vcf', help='Input .vcf file.')
    annparser.add_argument('ref_fasta', help='Reference .fasta file.')
    annparser.add_argument('bam', help='Input alignments.', action=CheckBam)
    annparser.add_argument('vcfout', help='Output .vcf.')
    annparser.add_argument('--chunk_size', default=500000, type=int,
        help='Chunk size for building pileups.')
    annparser.add_argument('--pad', default=25, type=int,
        help='Padding width either side of variant for realignment.')
    annparser.add_argument('--dpsp', default=False, action='store_true',
        help='Calulate depth and alignment score of spanning reads')

    # call variants by alignment of a consensus sequence to a reference
    c2vparser = toolsubparsers.add_parser('consensus2vcf',
        help='Call variants by alignment of a consensus sequence to a reference.',
        parents=[_log_level(), _align_chunking(), _region_ref_names()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    c2vparser.set_defaults(func=medaka.variant.vcf_from_fasta)
    c2vparser.add_argument('consensus', help='Input consensus .fasta file.')
    c2vparser.add_argument('ref_fasta', help='Reference .fasta file.')
    c2vparser.add_argument('--bam',
                           help='Existing bam file.')
    c2vparser.add_argument('--out_prefix', default='consensus2vcf',
                           help='Output prefix for .bam and .vcf files.')

    # convert a vcf to tsv
    tsvparser = toolsubparsers.add_parser('vcf2tsv',
        help='convert vcf to tsv, unpacking info and sample columns.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    tsvparser.set_defaults(func=medaka.vcf.vcf2tsv)
    tsvparser.add_argument('vcf', help='Input .vcf file.')

    # find homozygous and heterozygous regions in a VCF
    hzregparser = toolsubparsers.add_parser('homozygous_regions',
        help='Find homozygous regions from a diploid vcf.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    hzregparser.set_defaults(func=medaka.vcf.get_homozygous_regions)
    hzregparser.add_argument('vcf', help='Input .vcf file.')
    hzregparser.add_argument('region', help='Genomic region within which to find homozygous sub-regions.')
    hzregparser.add_argument('--min_len', type=int, default=1000,
                             help='Minimum region length.')
    hzregparser.add_argument('--suffix', help='Output suffix.', default='regions.txt')

    # resolve model and print full model file path to screen
    rmparser = toolsubparsers.add_parser('resolve_model',
        help='Resolve model and print full file path.',
        parents=[_log_level(), _model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rmparser.set_defaults(func=print_model_path)

    # check if feature encoder is RLE
    rleparser = toolsubparsers.add_parser('is_rle_model',
        help='Check if a model is an RLE model.',
        parents=[_model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rleparser.set_defaults(func=is_rle_model)

    # Request alignments parameters for model
    alignmentparser = toolsubparsers.add_parser('get_alignment_params',
        help='Get alignment parameters appropriate for a model',
        parents=[_model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    alignmentparser.set_defaults(func=print_alignment_params)

    # append RLE tags to a bam from hdf
    rlebamparser = toolsubparsers.add_parser('rlebam',
        description='Add RLE tags from HDF to bam. (input bam from stdin)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    rlebamparser.add_argument('read_index', help='Two column .tsv mapping read_ids to .hdf filepaths.')
    rlebamparser.add_argument('--workers', type=int, default=4, help='Number of worker processes.')
    rlebamparser.set_defaults(func=medaka.rle.rlebam)

    # print all model tags followed by default
    lmparser = toolsubparsers.add_parser('list_models',
        help='List all models.',
        parents=[_log_level(), _model_arg()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    lmparser.set_defaults(func=print_all_models)

    # write a bed file of regions spanned by a hdf feature / probs file.
    bdparser = toolsubparsers.add_parser('hdf_to_bed',
        help='Write a bed file of regions spanned by a hdf file.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    bdparser.set_defaults(func=medaka.variant.samples_to_bed)
    bdparser.add_argument('inputs', nargs='+',
                          help='Consensus .hdf files.')
    bdparser.add_argument('output', help='Output .bed.', default='medaka.bed')

    # resolve all available models, downloading missing ones.
    dwnldparser = toolsubparsers.add_parser('download_models',
        help='Download model files for any models not already installed.',
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    dwnldparser.add_argument(
        '--models', nargs='+', default=medaka.options.allowed_models,
        help='Model(s) to download to cache.')
    dwnldparser.set_defaults(func=download_models)

    return parser


def main():
    # Some users report setting this helps them resolve issues on their
    # filesystems.
    os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
    parser = medaka_parser()
    args = parser.parse_args()

    # https://github.com/tensorflow/tensorflow/issues/26691
    # we have local imports of tf generally, but this is pseudo global
    # however with tf1.14 tensorflow import time looks faster than previous
    # so we could/should stick to more standard imports
    import tensorflow as tf
    import absl.logging
    logging.root.removeHandler(absl.logging._absl_handler)
    absl.logging._warn_preinit_stderr = False

    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
    logger = logging.getLogger(__package__)
    logger.setLevel(args.log_level)

    if args.command == 'tools' and not hasattr(args, 'func'):
        # display help if given `medaka tools (--help)`
        # TODO: is there a cleaner way to access this?
        parser.__dict__['_actions'][1].choices['tools'].print_help()
    else:
        _validate_common_args(args)
        args.func(args)
