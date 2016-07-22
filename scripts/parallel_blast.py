#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from glob import glob
from os import makedirs, system
from os.path import exists, split, splitext, isfile
from subprocess import check_call, CalledProcessError

from bfillings.formatdb import build_blast_db_from_fasta_path

from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.util import load_qiime_config, get_options_lookup
from qiime.parallel.blast import ParallelBlaster

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel BLAST"""
script_info['script_description'] = """This script for performing blast while\
 making use of multicore/multiprocessor environments to perform analyses in\
 parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """BLAST $PWD/inseqs.fasta (-i) against a blast database created from\
 $PWD/refseqs.fasta (-r). Store the results in $PWD/blast_out/ (-o). ALWAYS\
 SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will\
 generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/inseqs.fasta -r $PWD/refseqs.fasta -o $PWD/blast_out/\
 -e 0.001"""))


script_info['output_description'] = """ """
script_info['required_options'] = [
    make_option('-i', '--infile_path', action='store',
                type='existing_filepath', dest='infile_path',
                help='Path of sequences to use as queries [REQUIRED]'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='name of output directory for blast jobs [REQUIRED]')
]
script_info['optional_options'] = [
    make_option('-c', '--disable_low_complexity_filter',
                default=False, action='store_true',
                help='disable filtering of low-complexity sequences '
                '(i.e., -F F is passed to blast) [default: %default]'),
    make_option('-e', '--e_value', action='store',
                type='float', default=1e-30, dest='e_value',
                help='E-value threshold for blasts [default: %default]'),
    make_option('-n', '--num_hits', action='store',
                type='int', default=1, dest='num_hits',
                help='number of hits per query for blast results [default: %default]'),
    make_option('-w', '--word_size', action='store',
                type='int', default=30, dest='word_size',
                help='word size for blast searches [default: %default]'),
    make_option('-a', '--blastmat_dir', action='store',
                type='string', help='full path to directory containing ' +
                'blastmat file [default: %default]',
                default=qiime_config['blastmat_dir']),
    make_option(
        '-r', '--refseqs_path', action='store', type='existing_filepath',
        help='Path to fasta sequences to search against. Required if ' +
             '-b is not provided.'),
    make_option('-b', '--blast_db', type='blast_db',
                help='Name of pre-formatted BLAST database. Required if ' +
                '-r is not provided.'),
    options_lookup['jobs_to_start'],
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if not (opts.refseqs_path or opts.blast_db):
        option_parser.error('Either a blast db (via -b) or a collection of '
                            'reference sequences (via -r) must be passed')
    if opts.refseqs_path and opts.blast_db:
        option_parser.error('You should provide only a blast db (via -b) '
                            'or a collection of reference sequences (via -r), but not both')

    # create dict of command-line options
    params = eval(str(opts))

    parallel_runner = ParallelBlaster(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)

    parallel_runner(opts.infile_path,
                    opts.output_dir,
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)


if __name__ == "__main__":
    main()
