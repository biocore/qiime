#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from glob import glob
from os.path import isfile
from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.beta_diversity import get_phylogenetic_metric
from qiime.beta_diversity import list_known_metrics
from qiime.util import load_qiime_config, get_options_lookup
from qiime.parallel.beta_diversity import (ParallelBetaDiversitySingle,
                                           ParallelBetaDiversityMultiple)

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel beta diversity"""
script_info['script_description'] = """This script performs like the\
 beta_diversity.py script, but is intended to make use of\
 multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Apply beta_diversity.py in parallel to multiple otu tables""",
     """Apply the unweighted_unifrac and weighted_unifrac metrics (modify with -m)\
 to all otu tables in rarefied_otu_tables (-i) and write the resulting output\
 files to bdiv/ (-o, will be created if it doesn't exist). Use the rep_set.tre\
 (-t) to compute phylogenetic diversity metrics. ALWAYS SPECIFY ABSOLUTE FILE\
 PATHS (absolute path represented here as $PWD, but will generally look\
 something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/rarefied_otu_tables/ -o $PWD/bdiv/ -t $PWD/rep_set.tre"""))

script_info['script_usage'].append(
    ("""Apply beta_diversity.py in parallel to a single otu table""", """ """,
     """%prog -i $PWD/otu_table.biom -o $PWD/bdiv_single/ -t $PWD/rep_set.tre"""))

script_info['output_description'] = """The output of %prog is a folder containing\
 text files, each a distance matrix between samples."""

script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_path',
                help='input path, must be directory [REQUIRED]'),
    make_option('-o', '--output_path', type='new_dirpath',
                help='output path, must be directory [REQUIRED]'),
]

script_info['optional_options'] = [
    make_option(
        '-m', '--metrics', default='unweighted_unifrac,weighted_unifrac',
        type='multiple_choice', mchoices=list_known_metrics(),
        help='Beta-diversity metric(s) to use. A comma-separated list should be' +
        ' provided when multiple metrics are specified. [default: %default]'),
    make_option('-t', '--tree_path', type='existing_filepath',
                help='path to newick tree file, required for phylogenetic metrics' +
                ' [default: %default]'),
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep'],
    options_lookup['jobs_to_start'],
    make_option('-f', '--full_tree', action="store_true",
                help='By default, each job removes calls _fast_unifrac_setup to remove\
 unused parts of the tree. pass -f if you already have a minimal tree, and\
 this script will run faster'),

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    params = eval(str(opts))

    params['metrics'] = ','.join(opts.metrics)

    # create local copies of command-line options
    input_path = opts.input_path
    output_dir = opts.output_path
    metrics_list = opts.metrics
    tree_fp = opts.tree_path

    # Check the tree exists if phylogenetically-aware measure is used
    for metric in metrics_list:
        try:
            metric_f = get_phylogenetic_metric(metric)
            if tree_fp is None:
                option_parser.error("metric %s requires a tree, but "
                                    "none found" % metric)
        except AttributeError:
            pass
    if isfile(input_path):
        # single otu table mode
        parallel_runner = ParallelBetaDiversitySingle(
            cluster_jobs_fp=opts.cluster_jobs_fp,
            jobs_to_start=opts.jobs_to_start,
            retain_temp_files=opts.retain_temp_files,
            suppress_polling=opts.suppress_polling,
            seconds_to_sleep=opts.seconds_to_sleep)
        parallel_runner(input_path,
                        output_dir,
                        params,
                        job_prefix=opts.job_prefix,
                        poll_directly=opts.poll_directly,
                        suppress_submit_jobs=opts.suppress_submit_jobs)

    else:
        input_fps = glob('%s/*' % input_path)
        parallel_runner = ParallelBetaDiversityMultiple(
            cluster_jobs_fp=opts.cluster_jobs_fp,
            jobs_to_start=opts.jobs_to_start,
            retain_temp_files=opts.retain_temp_files,
            suppress_polling=opts.suppress_polling,
            seconds_to_sleep=opts.seconds_to_sleep)
        parallel_runner(input_fps,
                        output_dir,
                        params,
                        job_prefix=opts.job_prefix,
                        poll_directly=opts.poll_directly,
                        suppress_submit_jobs=opts.suppress_submit_jobs)

if __name__ == "__main__":
    main()
