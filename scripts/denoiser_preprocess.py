#!/usr/bin/env python

"""Preprocess 454 sequencing data."""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
# remember to add yourself if you make changes
__credits__ = ["Jens Reeder", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from os.path import exists
from os import remove, rename, rmdir, makedirs
from random import sample
from tempfile import mkdtemp

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option
from qiime.denoiser.utils import files_exist, store_mapping
from qiime.denoiser.preprocess import preprocess

options_lookup = get_options_lookup()
STANDARD_BACTERIAL_PRIMER = "CATGCTGCCTCCCGTAGGAGT"

# denoiser_preprocess.py
script_info = {}
script_info[
    'brief_description'] = """Run phase of denoiser algorithm: prefix clustering"""
script_info['script_description'] = """The script denoiser_preprocess.py runs the first clustering phase
which groups reads based on common prefixes."""

script_info['script_usage'] = [
    ("",
     """Run program on flowgrams in 454Reads.sff. Remove reads which are not in split_lib_filtered_seqs.fasta.
Remove primer CATGCTGCCTCCCGTAGGAGT from reads before running phase I""",
     """%prog -i Fasting_Example.sff.txt -f seqs.fna -p CATGCTGCCTCCCGTAGGAGT """ )
]

script_info['output_description'] = """
prefix_dereplicated.sff.txt: human readable sff file containing the flowgram of the
                             cluster representative of each cluster.

prefix_dereplicated.fasta: Fasta file containing the cluster representative of each cluster.

prefix_mapping.txt: This file contains the actual clusters. The cluster centroid is given first,
                    the cluster members follw after the ':'.
"""

script_info['required_options'] = [

    make_option('-i', '--input_files', action='store',
                type='existing_filepaths', dest='sff_fps',
                help='path to flowgram files (.sff.txt), ' +
                'comma separated')
]

script_info['optional_options'] = [
    make_option('-f', '--fasta_file', action='store', type='string',
                dest='fasta_fp', help='path to fasta input file ' +
                '[default: %default]', default=None),

    make_option('-s', '--squeeze', action='store_true', dest='squeeze',
                help='Use run-length encoding for prefix ' +
                'filtering [default: %default]', default=False),

    make_option('-l', '--log_file', action='store',
                type='string', dest='log_fp', help='path to log file ' +
                '[default: %default]', default="preprocess.log"),

    make_option('-p', '--primer', action='store',
                type='string', dest='primer', help='primer sequence ' +
                'used for the amplification [default: %default]',
                default=STANDARD_BACTERIAL_PRIMER),

    make_option('-o', '--output_dir', action='store',
                type='string', dest='output_dir',
                help='path to output directory ' +
                '[default: %default]', default="/tmp/")
]

script_info['version'] = __version__


def main(commandline_args=None):
    parser, opts, args = parse_command_line_parameters(**script_info)

    # make tmp and output dir
    try:
        tmp_dir = mkdtemp(dir=opts.output_dir, suffix="/")
    except OSError:
        exit("Creating temporary directory failed")
    if(not exists(opts.output_dir)):
        try:
            makedirs(opts.output_dir)
        except OSError:
            exit("Creating output directory failed")

    # open logger
    log_fh = None
    if opts.verbose:
        # append to the log file of the master process
        log_fh = open(opts.output_dir + "/" + opts.log_fp, "a", 0)
        log_fh.write("SFF files: %s" % ', '.join(opts.sff_fps))
        log_fh.write("Fasta file: %s\n" % opts.fasta_fp)
        log_fh.write("Output dir: %s\n" % opts.output_dir)
        log_fh.write("Squeeze Seqs: %s\n" % opts.squeeze)
        log_fh.write("Primer sequence: %s\n" % opts.primer)

    (deprefixed_sff_fp, l, mapping, seqs) = \
        preprocess(opts.sff_fps, log_fh, fasta_fp=opts.fasta_fp,
                   out_fp=tmp_dir,
                   verbose=opts.verbose, squeeze=opts.squeeze,
                   primer=opts.primer)

    # explicitly close log file, as this file can be shared with the master
    # Closing it here assures that all preprocess writes happen before the
    # master writes
    if log_fh:
        log_fh.close()

    # move files to output dir
    rename(tmp_dir + "/prefix_dereplicated.sff.txt",
           opts.output_dir + "/prefix_dereplicated.sff.txt")
    rename(tmp_dir + "/prefix_dereplicated.fasta",
           opts.output_dir + "/prefix_dereplicated.fasta")
    rename(
        tmp_dir +
        "/prefix_mapping.txt",
        opts.output_dir +
        "/prefix_mapping.txt")
    rmdir(tmp_dir)

if __name__ == "__main__":
    main()
