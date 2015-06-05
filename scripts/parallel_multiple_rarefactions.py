#!/usr/bin/env python
# File created on 14 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import (parse_command_line_parameters, make_option)
from os.path import split, splitext, join
from qiime.util import get_options_lookup
from qiime.parallel.multiple_rarefactions import ParallelMultipleRarefactions

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel multiple file rarefaction"""
script_info[
    'script_description'] = """This script performs like the multiple_rarefactions.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""OTU tables of different depths""",
     """Build rarefied otu tables containing 10 (-m) to 140 (-x) sequences in steps of 10 (-s) with 2 (-n) repetions per number of sequences, from otu_table.biom (-i). Write the output files to the rarefied_otu_tables directory (-o, will be created if it doesn't exist). The name of the output files will be of the form rarefaction_<num_seqs>_<reptition_number>.biom. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -o $PWD/rarefied_otu_tables/ -m 10 -x 140 -s 10 -n 2 -i $PWD/otu_table.biom"""))

script_info['script_usage'].append(
    ("""OTU tables of the same depth""",
     """Build 8 rarefied otu tables each containing exactly 100 sequences per sample (even depth rarefaction). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -o $PWD/even_otu_tables/ -m 100 -x 100 -n 8 -i $PWD/otu_table.biom"""))

script_info[
    'output_description'] = """The result of parallel_multiple_rarefactions.py consists of a number of files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.biom, and are named in the following way: rarefaction_100_0.txt, where "100" corresponds to the sequences per sample and "0" for the iteration."""

script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_filepath',
                help='input filepath, (the otu table) [REQUIRED]'),
    make_option('-o', '--output_path', type='new_dirpath',
                help="write output rarefied otu tables here makes dir if it doesn't exist [REQUIRED]"),
    make_option('-m', '--min', type=int, help='min seqs/sample [REQUIRED]'),
    make_option('-x', '--max', type=int,
                help='max seqs/sample (inclusive) [REQUIRED]'),

]
script_info['optional_options'] = [
    make_option('-n', '--num_reps', dest='num_reps', default=10, type=int,
                help='num iterations at each seqs/sample level [default: %default]'),
    make_option(
        '--suppress_lineages_included', default=False, action="store_true",
        help='Exclude taxonomic (lineage) information for each OTU.'),
    make_option('-s', '--step', type=int, default=1,
                help='levels: min, min+step... for level <= max [default: %default]'),
    make_option('--subsample_multinomial', default=False, action='store_true',
                help='subsample using subsampling with replacement [default: %default]'),
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep'],
    options_lookup['jobs_to_start']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # create dict of command-line options
    params = eval(str(opts))

    if not opts.step > 0:
        option_parser.error(("Error: step size must be greater than 0.\n"
                             "If min = max, just leave step size at 1."))

    parallel_runner = ParallelMultipleRarefactions(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)
    parallel_runner(opts.input_path,
                    opts.output_path,
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)


if __name__ == "__main__":
    main()
