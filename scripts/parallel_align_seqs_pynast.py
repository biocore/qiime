#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup, make_option,
                        load_qiime_config)
from qiime.align_seqs import pairwise_alignment_methods
from qiime.parallel.align_seqs import ParallelAlignSeqsPyNast

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Parallel sequence alignment using PyNAST"""
script_info[
    'script_description'] = """A wrapper for the align_seqs.py PyNAST option, intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Align the input file (-i) against using PyNAST and write the output (-o) to $PWD/pynast_aligned_seqs/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/inseqs.fasta -o $PWD/pynast_aligned_seqs/"""))
script_info[
    'output_description'] = """This results in a multiple sequence alignment (FASTA-formatted)."""

script_info['required_options'] = [
    options_lookup['fasta_as_primary_input'],
    options_lookup['output_dir']
]

pairwise_alignment_method_choices = pairwise_alignment_methods.keys()
blast_db_default_help =\
    qiime_config['pynast_template_alignment_blastdb'] or \
    'created on-the-fly from template_alignment'

script_info['optional_options'] = [
    make_option('-t', '--template_fp', type='existing_filepath',
                help='Filepath for template alignment [default: %default]',
                default=qiime_config['pynast_template_alignment_fp']),
    make_option('-a', '--pairwise_alignment_method',
                type='choice', help='Method to use for pairwise alignments' +
                ' [default: %default]',
                default='uclust', choices=pairwise_alignment_method_choices),
    make_option('-d', '--blast_db', type='blast_db',
                dest='blast_db', help='Database to blast against' +
                ' [default: %s]' % blast_db_default_help,
                default=qiime_config['pynast_template_alignment_blastdb']),
    make_option('-e', '--min_length',
                type='int', help='Minimum sequence ' +
                'length to include in alignment [default: 75% of the' +
                ' median input sequence length]',
                default=-1),
    make_option('-p', '--min_percent_id', action='store',
                type='float', help='Minimum percent ' +
                'sequence identity to closest blast hit to include sequence in' +
                ' alignment, expressed as a real number between 0 and 100 '
                '[default: %default]', default=75.0),
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

    # create dict of command-line options
    params = eval(str(opts))

    parallel_runner = ParallelAlignSeqsPyNast(
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
