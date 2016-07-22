#!/usr/bin/env python
# File created on 09 Oct 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import getenv
from os.path import join
from qiime.util import (get_options_lookup, load_qiime_config, make_option,
                         parse_command_line_parameters)
from qiime.parallel.assign_taxonomy import ParallelUclustConsensusTaxonomyAssigner

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Parallel taxonomy assignment using the uclust consensus taxonomy assignment"""
script_info[
    'script_description'] = """This script performs like the assign_taxonomy.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Assign taxonomy to all sequences in the input file (-i) using the uclust consensus taxonomy assigner and write the results (-o) to $PWD/uclust_assigned_taxonomy/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/inseqs.fasta -o $PWD/uclust_assigned_taxonomy/"""))
script_info[
    'output_description'] = """Mapping of sequence identifiers to taxonomy and quality information."""
script_info['required_options'] = [
    make_option('-i', '--input_fasta_fp', action='store',
                type='existing_filepath',
                help='full path to fasta file containing query sequences [REQUIRED]'),
    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', help='path to store output files [REQUIRED]'),
]

default_reference_seqs_fp = qiime_config['assign_taxonomy_reference_seqs_fp']
default_id_to_taxonomy_fp = qiime_config['assign_taxonomy_id_to_taxonomy_fp']

script_info['optional_options'] = [
    make_option('-t', '--id_to_taxonomy_fp', action='store',
                type='existing_filepath', help='full path to '
                'id_to_taxonomy mapping file [default: %default]',
                default=default_id_to_taxonomy_fp),
    make_option('-r', '--reference_seqs_fp', action='store',
                help='Ref seqs to search against. [default: %default]',
                default=default_reference_seqs_fp, type='existing_filepath'),
    make_option('--min_consensus_fraction', type='float',
                help=('Minimum fraction of database hits that must have a '
                      'specific taxonomic assignment to assign that taxonomy '
                      'to a query [default: %default]'), default=0.51),
    make_option('--similarity', type='float',
                help=('Minimum percent similarity to consider a database '
                      'match a hit, expressed as a fraction between 0 and 1 '
                      '[default: %default]'), default=0.90),
    make_option('--uclust_max_accepts', type='int',
                help=('Number of database hits to consider when making '
                      'an assignment [default: %default]'), default=3),
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

    parallel_runner = ParallelUclustConsensusTaxonomyAssigner(
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
