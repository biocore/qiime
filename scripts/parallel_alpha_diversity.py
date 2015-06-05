#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import parse_command_line_parameters, make_option
from glob import glob
from os.path import join
from qiime.util import get_options_lookup
from qiime.parallel.alpha_diversity import ParallelAlphaDiversity
from qiime.alpha_diversity import list_known_metrics

script_info = {}
script_info['brief_description'] = """Parallel alpha diversity"""
script_info['script_description'] = """This script performs like the\
 alpha_diversity.py script, but is intended to make use of\
 multicore/multiprocessor environments to perform analyses in parallel."""

script_info['script_usage'] = []

script_info['script_usage'].append(("""Example""",
                                    """Apply the observed_OTUs, chao1, PD_whole_tree metrics (-m) to all otu\
 tables in rarefied_otu_tables/ (-i) and write the resulting output files to\
 adiv/ (-o, will be created if it doesn't exist). Use the rep_set.tre (-t) to\
 compute phylogenetic diversity metrics. ALWAYS SPECIFY ABSOLUTE FILE PATHS\
 (absolute path represented here as $PWD, but will generally look something\
 like /home/ubuntu/my_analysis/).""",
                                    """%prog -i $PWD/rarefied_otu_tables -o $PWD/adiv\
 -m observed_otus,chao1,PD_whole_tree -t $PWD/rep_set.tre"""))

script_info['output_description'] = """The resulting output will be the same\
 number of files as supplied by the user. The resulting files are tab-delimited\
 text files, where the columns correspond to alpha diversity metrics and the\
 rows correspond to samples and their calculated diversity measurements. """

script_info['version'] = __version__

options_lookup = get_options_lookup()

script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_dirpath',
                help='input path, must be directory [REQUIRED]'),
    make_option('-o', '--output_path', type='new_dirpath',
                help='output path, must be directory [REQUIRED]'),
]

script_info['optional_options'] = [
    make_option('-t', '--tree_path', type='existing_filepath',
                help='path to newick tree file, required for phylogenetic metrics' +
                ' [default: %default]'),
    make_option('-m', '--metrics', type='multiple_choice',
                mchoices=list_known_metrics(),
                help='metrics to use, comma delimited',
                default='PD_whole_tree,chao1,observed_otus'),
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep'],
    options_lookup['jobs_to_start']
]


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    params = eval(str(opts))

    params['metrics'] = ','.join(opts.metrics)

    parallel_runner = ParallelAlphaDiversity(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)
    input_fps = glob(join(opts.input_path, '*'))
    parallel_runner(input_fps,
                    opts.output_path,
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)

if __name__ == "__main__":
    main()
