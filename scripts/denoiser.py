#!/usr/bin/env python
from __future__ import division

""" A routine to clean up 454 sequencing data."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os import makedirs, remove
from os.path import exists
from tempfile import mkdtemp

from skbio.util import create_dir
from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option

from qiime.denoiser.preprocess import STANDARD_BACTERIAL_PRIMER
from qiime.denoiser.flowgram_clustering import denoise_seqs, denoise_per_sample
from qiime.denoiser.utils import files_exist, get_denoiser_data_dir,\
    cat_sff_files

options_lookup = get_options_lookup()
DENOISER_DATA_DIR = get_denoiser_data_dir()
RELATIVE_DENOISER_DATA_DIR = '<qiime-install-path>/qiime/support_files/denoiser/Data/'

# denoiser.py
script_info = {}
script_info['brief_description'] = """Remove noise from  454 sequencing data"""
script_info[
    'script_description'] = """The denoiser removes sequencing noise characteristic to pyrosequencing by flowgram clustering. For a detailed explanation of the underlying algorithm see (Reeder and Knight, Nature Methods 7(9), 2010)."""

script_info['script_usage'] = [
    ("",
     """Run denoiser on flowgrams in 454Reads.sff.txt with read-to-barcode mapping in seqs.fna,
put results into Outdir, log progress in Outdir/denoiser.log""",
     """%prog -i 454Reads.sff.txt -f seqs.fna -v -o Outdir"""),

    ("Multiple sff.txt files",
     """Run denoiser on two flowgram files in 454Reads_1.sff.txt and 454Reads_2.sff.txt
with read-to-barcode mapping in seqs.fna, put results into Outdir,
log progress in Outdir/denoiser.log""",
     """%prog -i 454Reads_1.sff.txt,454Reads_2.sff.txt -f seqs.fna -v -o Outdir"""),

    ("Denoise multiple library separately",
     """Run denoiser on flowgrams in 454Reads.sff.txt with read-to-barcode mapping in seqs.fna,
split input files into libraries and process each library separately,
put results into Outdir, log progress in Outdir/denoiser.log""",
     "%prog -S -i 454Reads.sff.txt -f seqs.fna -v -o Outdir"),

    ("Resuming a failed run",
     """Resume a previous denoiser run from breakpoint stored in Outdir_from_failed_run/checkpoints/checkpoint100.pickle.
The checkpoint option requires the -p or --preprocess option, which usually can be set to the output dir of the failed run.
All other arguments must be identical to the failed run.""",
     "%prog -i 454Reads.sff.txt -f seqs.fna -v -o Outdir_resumed -p Outdir_from_failed_run --checkpoint Outdir_from_failed_run/checkpoints/checkpoint100.pickle"),
]

script_info['output_description'] = """

centroids.fasta: The cluster representatives of each cluster

singletons.fasta: contains all unclustered reads

denoiser_mapping.txt: This file contains the actual clusters. The cluster centroid is given first,
                    the cluster members follow after the ':'.

checkpoints/ : directory with checkpoints

Note that the centroids and singleton files are disjoint. For most downstream analyses one wants to cat the two files.
"""

script_info['required_options'] = [

    make_option('-i', '--input_files', action='store',
                type='existing_filepaths', dest='sff_fps',
                help='path to flowgram files (.sff.txt), ' +
                'comma separated')
]

script_info['optional_options'] = [

    make_option('-f', '--fasta_fp', action='store',
                type='string', dest='fasta_fp',
                help='path to fasta input file. ' +
                'Reads not in the fasta file are filtered out ' +
                'before denoising. File format is as produced by ' +
                'split_libraries.py ' +
                '[default: %default]',
                default=None),

    make_option('-o', '--output_dir', action='store',
                type='string', dest='output_dir', help='path to output' +
                ' directory [default: random dir in ./]',
                default=None),

    make_option('-c', '--cluster', action='store_true',
                dest='cluster', help='Use cluster/multiple CPUs for ' +
                'flowgram alignments [default: %default]',
                default=False),

    make_option('-p', '--preprocess_fp', action='store',
                type='string', dest='preprocess_fp',
                help='Do not do preprocessing (phase I),'
                + 'instead use already preprocessed data in PREPROCESS_FP',
                default=None),

    make_option('--checkpoint_fp', action='store',
                type='string', dest='checkpoint_fp',
                help='Resume denoising from checkpoint. '
                + 'Be careful when changing parameters for'
                + ' a resumed run. Requires -p option. ' +
                ' [default: %default]',
                default=None),

    make_option('-s', '--squeeze', action='store_true',
                dest='squeeze', help='Use run-length encoding for prefix ' +
                'filtering in phase I [default: %default]',
                default=False),

    make_option('-S', '--split', action='store_true',
                dest='split', help='Split input into per library sets ' +
                'and denoise separately [default: %default]',
                default=False),

    make_option('--force', action='store_true',
                dest='force', help='Force overwrite of existing '
                + 'directory [default: %default]',
                default=False),

    make_option('--primer', action='store',
                type='string', dest='primer',
                help='primer sequence ' +
                '[default: %default]',
                default=STANDARD_BACTERIAL_PRIMER),

    make_option('-n', '--num_cpus', action='store',
                type='int', dest='num_cpus',
                help='number of cpus, requires -c ' +
                '[default: %default]', default=1),

    make_option('-m', '--max_num_iterations', action='store',
                type='int', dest='max_num_iter',
                help='maximal number of iterations in phase II. ' +
                'None means unlimited iterations ' +
                '[default: %default]', default=None),

    make_option('-b', '--bail_out', action='store',
                type='int', dest='bail',
                help='stop clustering in phase II with ' +
                'clusters smaller or equal than BAILde' +
                ' [default: %default]', default=1),

    make_option('--percent_id', action='store',
                type='float', dest='percent_id',
                help='sequence similarity clustering '
                'threshold, expressed as a fraction between 0 and 1 '
                '[default: %default]', default=0.97),

    make_option('--low_cut_off', action='store',
                type='float', dest='low_cutoff',
                help='low clustering threshold for phase II ' +
                '[default: %default]', default=3.75),

    make_option('--high_cut_off', action='store',
                type='float', dest='high_cutoff',
                help='high clustering threshold for phase III ' +
                '[default: %default]', default=4.5),

    make_option('--low_memory', action='store_true',
                dest='low_memory', help='Use slower, low ' +
                'memory method [default: %default]',
                default=False),

    make_option('-e', '--error_profile', action='store',
                dest='error_profile', help='path to error profile ' +
                '[default= %s]' % (RELATIVE_DENOISER_DATA_DIR +
                                   'FLX_error_profile.dat'),
                default=DENOISER_DATA_DIR + 'FLX_error_profile.dat'),

    # might be needed once we switch to Titanium as default
    #    make_option('-flx', action='store_true',
    #                      dest='flx', help='shortcut for '+
    #                      "-e %s/FLX_error_profile.dat --low_cut_off=3.75 --high_cut_off=4.5" % DENOISER_DATA_DIR),

    make_option('--titanium', action='store_true',
                dest='titanium', help='shortcut for '
                '-e ' + RELATIVE_DENOISER_DATA_DIR +
                '/Titanium_error_profile.dat --low_cut_off=4 --high_cut_off=5 . ' +
                'Warning: overwrites all previous cut-off values ' +
                '[DEFAULT: %default]', default=False)
]

script_info['version'] = __version__


def main(commandline_args=None):
    parser, opts, args = parse_command_line_parameters(**script_info)

    if(opts.checkpoint_fp):
        bp_fp = opts.checkpoint_fp
        if not exists(bp_fp):
            parser.error(
                'Specified checkpoint file does not exist: %s' %
                bp_fp)

    # peek into sff.txt files to make sure they are parseable
    # cat_sff_fles is lazy and only reads header
    flowgrams, header = cat_sff_files(map(open, opts.sff_fps))

    if(opts.split and opts.preprocess_fp):
        parser.error('Options --split and --preprocess_fp are exclusive')

    if(opts.preprocess_fp):
        pp_fp = opts.preprocess_fp
        if not exists(opts.preprocess_fp):
            parser.error(
                'Specified preprocess directory does not exist: %s' %
                opts.preprocess_fp)
        if not files_exist('%s/prefix_mapping.txt,%s/prefix_dereplicated.fasta' % (pp_fp, pp_fp)):
            parser.error('Specified preprocess directory does not contain expected files: ' +
                         'prefix_mapping.txt and prefix_dereplicated.fasta')

    if opts.titanium:
        opts.error_profile = DENOISER_DATA_DIR + 'Titanium_error_profile.dat'
        opts.low_cutoff = 4
        opts.high_cutoff = 5

    if not exists(opts.error_profile):
        parser.error(
            'Specified error profile %s does not exist' %
            opts.error_profile)

    if opts.output_dir:
        # make sure it always ends on /
        tmpoutdir = opts.output_dir + "/"
        create_dir(tmpoutdir, not opts.force)
    else:
        # make random dir in current dir
        tmpoutdir = mkdtemp(dir="", prefix="denoiser_", suffix="/")


    log_fp = 'denoiser.log'

    if opts.split:
        denoise_per_sample(
            opts.sff_fps, opts.fasta_fp, tmpoutdir, opts.cluster,
            opts.num_cpus, opts.squeeze, opts.percent_id, opts.bail,
            opts.primer, opts.low_cutoff, opts.high_cutoff, log_fp,
            opts.low_memory, opts.verbose, opts.error_profile, opts.max_num_iter,
            opts.titanium)
    else:
        denoise_seqs(
            opts.sff_fps, opts.fasta_fp, tmpoutdir, opts.preprocess_fp, opts.cluster,
            opts.num_cpus, opts.squeeze, opts.percent_id, opts.bail, opts.primer,
            opts.low_cutoff, opts.high_cutoff, log_fp, opts.low_memory,
            opts.verbose, opts.error_profile, opts.max_num_iter, opts.titanium,
            opts.checkpoint_fp)

if __name__ == "__main__":
    main()
