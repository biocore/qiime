#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Dan Knights", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from qiime.parallel.pick_otus import ParallelPickOtusBlast

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel pick otus using BLAST"""
script_info[
    'script_description'] = """This script performs like the pick_otus.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Pick OTUs by blasting $PWD/inseqs.fasta against $PWD/refseqs.fasta and write the output to the $PWD/blast_otus/ directory. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/blast_otus/"""))
script_info[
    'output_description'] = """The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file."""

script_info['required_options'] = [
    make_option('-i', '--input_fasta_fp', action='store',
                type='existing_filepath', help='full path to ' +
                'input_fasta_fp'),
    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', help='path to store output files')
]

script_info['optional_options'] = [
    make_option('-e', '--max_e_value',
                help='Max E-value ' +
                '[default: %default]', default='1e-10'),

    make_option('-s', '--similarity', action='store',
                type='float', help='Sequence similarity ' +
                'threshold [default: %default]', default=0.97),

    make_option('-r', '--refseqs_fp', action='store',
                type='existing_filepath', help='full path to ' +
                'template alignment [default: %default]'),

    make_option('-b', '--blast_db', action='store',
                type='blast_db', help='database to blast against ' +
                '[default: %default]'),

    make_option('--min_aligned_percent',
                help=('Minimum percent of query sequence that can be aligned '
                      'to consider a hit, expressed as a fraction between 0 '
                      'and 1 (BLAST OTU picker only) [default: %default]'),
                default=0.50, type='float'),

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

    if opts.blast_db is None and opts.refseqs_fp is None:
        option_parser.error('Either blast_db or refseqs_fp must be provided.')

    # create dict of command-line options
    params = eval(str(opts))

    parallel_runner = ParallelPickOtusBlast(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)
    parallel_runner(opts.input_fasta_fp,
                    opts.output_dir,
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)


if __name__ == "__main__":
    main()
