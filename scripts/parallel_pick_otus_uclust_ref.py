#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from qiime.parallel.pick_otus import ParallelPickOtusUclustRef


############################
# Script functionality

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel pick otus using uclust_ref"""
script_info[
    'script_description'] = """This script works like the pick_otus.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Pick OTUs by searching $PWD/inseqs.fasta against $PWD/refseqs.fasta with reference-based uclust and write the output to the $PWD/blast_otus/ directory. This is a closed-reference OTU picking process. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/ucref_otus/"""))
script_info[
    'output_description'] = """The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an OTU identifier which is the reference sequence identifier, and the remaining fields correspond to sequence identifiers assigned to that OTU. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file."""

script_info['required_options'] = [
    make_option('-i', '--input_fasta_fp', action='store',
                type='existing_filepath', help='full path to ' +
                'input_fasta_fp'),

    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', help='path to store output files'),

    make_option('-r', '--refseqs_fp', action='store',
                type='existing_filepath', help='full path to ' +
                'reference collection')
]

script_info['optional_options'] = [

    make_option('-s', '--similarity', action='store',
                type='float', help='Sequence similarity ' +
                'threshold [default: %default]', default=0.97),
    make_option('-z', '--enable_rev_strand_match', action='store_true',
                default=False,
                help=('Enable reverse strand matching for uclust otu picking, '
                      'will double the amount of memory used. [default: %default]')),
    make_option('-A', '--optimal_uclust', action='store_true',
                default=False,
                help=('Pass the --optimal flag to uclust for uclust otu'
                      ' picking. [default: %default]')),
    make_option('-E', '--exact_uclust', action='store_true',
                default=False,
                help=('Pass the --exact flag to uclust for uclust otu'
                      ' picking. [default: %default]')),
    make_option('--max_accepts', type='int', default=1,
                help="max_accepts value to uclust and uclust_ref [default: %default]"),
    make_option('--max_rejects', type='int', default=8,
                help="max_rejects value to uclust and uclust_ref [default: %default]"),
    make_option('--stepwords', type='int', default=8,
                help="stepwords value to uclust and uclust_ref [default: %default]"),
    make_option('--word_length', type='int', default=8,
                help="w value to uclust and uclust_ref [default: %default]"),
    make_option('--uclust_stable_sort', default=True, action='store_true',
                help=("Deprecated: stable sort enabled by default, pass "
                      "--uclust_suppress_stable_sort to disable [default: %default]")),
    make_option(
        '--suppress_uclust_stable_sort', default=False, action='store_true',
        help=(
            "Don't pass --stable-sort to uclust [default: %default]")),
    make_option('-d', '--save_uc_files', default=True, action='store_false',
                help=("Enable preservation of intermediate uclust (.uc) files "
                      "that are used to generate clusters via uclust. "
                      "[default: %default]")),
    make_option('--denovo_otu_id_prefix', default=None,
                help=("OTU identifier prefix (string) for the de novo uclust"
                      " OTU picker [default: %default, OTU ids are ascending"
                      " integers]")),

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
    params['stable_sort'] = not opts.suppress_uclust_stable_sort

    parallel_runner = ParallelPickOtusUclustRef(
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
