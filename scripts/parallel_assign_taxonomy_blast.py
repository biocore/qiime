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

from os.path import join, abspath
from qiime.util import (get_options_lookup, load_qiime_config, make_option,
                         parse_command_line_parameters)
from qiime.parallel.assign_taxonomy import ParallelBlastTaxonomyAssigner

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

default_reference_seqs_fp = qiime_config['assign_taxonomy_reference_seqs_fp']
default_id_to_taxonomy_fp = qiime_config['assign_taxonomy_id_to_taxonomy_fp']

script_info = {}

script_info[
    'brief_description'] = """Parallel taxonomy assignment using BLAST"""
script_info[
    'script_description'] = """This script performs like the assign_taxonomy.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""

script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Assign taxonomy to all sequences in the input file (-i) using BLAST with the id to taxonomy mapping file (-t) and reference sequences file (-r), and write the results (-o) to $PWD/blast_assigned_taxonomy/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/inseqs.fasta -t $PWD/id_to_tax.txt -r $PWD/refseqs.fasta -o $PWD/blast_assigned_taxonomy/"""))

script_info[
    'output_description'] = """Mapping of sequence identifiers to taxonomy and quality scores."""

script_info['required_options'] = [
    make_option('-i', '--input_fasta_fp',
                type='existing_filepath', help='full path to ' +
                'input_fasta_fp [REQUIRED]'),
    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', help='full path to store output files ' +
                '[REQUIRED]')
]

script_info['optional_options'] = [
    make_option('-r', '--reference_seqs_fp', type='existing_filepath',
                help='Ref seqs to blast against.  Must provide either --blast_db or '
                '--reference_seqs_db for assignment with blast [default: %s]'
                % default_reference_seqs_fp,
                default=default_reference_seqs_fp),
    make_option('-b', '--blast_db', type='blast_db',
                help='Database to blast against.  Must provide either --blast_db or '
                '--reference_seqs_db for assignment with blast [default: %default]'),
    make_option('-e', '--e_value', type='float',
                help='Maximum e-value to record an assignment, only used for blast '
                'method [default: %default]', default=0.001),
    make_option('-B', '--blastmat_dir', action='store',
                type='string', help='full path to directory containing ' +
                'blastmat file [default: %default]',
                default=qiime_config['blastmat_dir']),
    options_lookup['jobs_to_start'],
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep']
]

if default_id_to_taxonomy_fp:
    script_info['optional_options'].append(
        make_option('-t', '--id_to_taxonomy_fp', action='store',
                    type='existing_filepath', help='full path to ' +
                    'id_to_taxonomy mapping file [default: %s]' % default_id_to_taxonomy_fp,
                    default=default_id_to_taxonomy_fp))
else:
    script_info['required_options'].append(
        make_option('-t', '--id_to_taxonomy_fp', action='store',
                    type='existing_filepath', help='full path to ' +
                    'id_to_taxonomy mapping file [REQUIRED]'))

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if not (opts.reference_seqs_fp or opts.blast_db):
        option_parser.error('Either a blast db (via -b) or a collection of '
                            'reference sequences (via -r) must be passed to '
                            'assign taxonomy using blast.')

    # create dict of command-line options
    params = eval(str(opts))

    parallel_runner = ParallelBlastTaxonomyAssigner(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)

    parallel_runner(opts.input_fasta_fp,
                    abspath(opts.output_dir),
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)


if __name__ == "__main__":
    main()
