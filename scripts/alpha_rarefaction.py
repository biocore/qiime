#!/usr/bin/env python
# File created on 04 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from os import makedirs
from qiime.util import load_qiime_config
from qiime.parse import parse_qiime_parameters
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 print_to_stdout,
                                 no_status_updates,
                                 validate_and_set_jobs_to_start)
from qiime.workflow.downstream import run_alpha_rarefaction

qiime_config = load_qiime_config()

options_lookup = get_options_lookup()
script_info = {}
script_info[
    'brief_description'] = """A workflow script for performing alpha rarefaction"""
script_info['script_description'] = """
The steps performed by this script are: Generate rarefied OTU tables; compute alpha diversity metrics for each rarefied OTU table; collate alpha diversity results; and generate alpha rarefaction plots."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Given an OTU table, a phylogenetic tree, a mapping file, and a max sample depth, compute alpha rarefaction plots for the PD, observed species and chao1 metrics. To specify alternative metrics pass a parameter file via -p. We generally recommend that the max depth specified here (-e) is the same as the even sampling depth provided to beta_diversity_through_plots (also -e). """,
     """%prog -i otu_table.biom -o arare_max100/ -t rep_set.tre -m Fasting_Map.txt -e 100"""))

script_info[
    'output_description'] = """The primary interface for the results will be OUTPUT_DIR/alpha_rarefaction_plots/rarefaction_plots.html where OUTPUT_DIR is the value you specify with -o.  You can open this in a web browser for interactive alpha rarefaction plots."""

script_info['script_usage_output_to_remove'] = ['arare_max100']

script_info['required_options'] = [
    make_option('-i', '--otu_table_fp', type='existing_filepath',
                help='the input otu table [REQUIRED]'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='path to the mapping file [REQUIRED]'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output directory [REQUIRED]'),
]
script_info['optional_options'] = [
    make_option('-p', '--parameter_fp', type='existing_filepath',
                help='path to the parameter file, which specifies changes' +
                ' to the default behavior. ' +
                'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .' +
                ' [if omitted, default values will be used]'),
    make_option('-n', '--num_steps', type='int',
                help='number of steps (or rarefied OTU table sizes) to make between ' +
                'min and max counts [default: %default]', default=10),
    make_option('-f', '--force', action='store_true',
                dest='force', help='Force overwrite of existing output directory' +
                ' (note: existing files in output_dir will not be removed)' +
                ' [default: %default]'),
    make_option('-w', '--print_only', action='store_true',
                dest='print_only', help='Print the commands but don\'t call them -- ' +
                'useful for debugging [default: %default]', default=False),
    make_option('-a', '--parallel', action='store_true',
                dest='parallel', default=False,
                help='Run in parallel where available [default: %default]'),
    make_option('-t', '--tree_fp', type='existing_filepath',
                help='path to the tree file [default: %default; ' +
                'REQUIRED for phylogenetic measures]'),
    make_option('--min_rare_depth', type='int',
                help='the lower limit of rarefaction depths ' +
                '[default: %default]', default=10),
    make_option('-e', '--max_rare_depth', type='int',
                help='the upper limit of rarefaction depths ' +
                '[default: median sequence/sample count]'),
    options_lookup['jobs_to_start_workflow'],
    make_option('--retain_intermediate_files', action='store_true', help='retain '
                'intermediate files: rarefied OTU tables (rarefaction) and alpha diversity '
                'results (alpha_div). By default these will be erased [default: %default]',
                default=False),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    verbose = opts.verbose

    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    num_steps = opts.num_steps
    verbose = opts.verbose
    print_only = opts.print_only
    parallel = opts.parallel
    min_rare_depth = opts.min_rare_depth
    max_rare_depth = opts.max_rare_depth
    retain_intermediate_files = opts.retain_intermediate_files

    if opts.parameter_fp:
        try:
            parameter_f = open(opts.parameter_fp, 'U')
        except IOError:
            raise IOError("Can't open parameters file (%s). Does it exist? Do you have read access?"
                          % opts.parameter_fp)
        params = parse_qiime_parameters(parameter_f)
        parameter_f.close()
    else:
        params = parse_qiime_parameters([])
        # empty list returns empty defaultdict for now

    jobs_to_start = opts.jobs_to_start
    default_jobs_to_start = qiime_config['jobs_to_start']
    validate_and_set_jobs_to_start(params,
                                   jobs_to_start,
                                   default_jobs_to_start,
                                   parallel,
                                   option_parser)

    try:
        makedirs(output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # Since the analysis can take quite a while, I put this check
            # in to help users avoid overwriting previous output.
            option_parser.error("Output directory already exists. Please choose"
                                " a different directory, or force overwrite with -f.")

    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially

    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    run_alpha_rarefaction(otu_table_fp=otu_table_fp,
                          mapping_fp=mapping_fp,
                          output_dir=output_dir,
                          command_handler=command_handler,
                          params=params,
                          qiime_config=qiime_config,
                          tree_fp=tree_fp,
                          num_steps=num_steps,
                          parallel=parallel,
                          min_rare_depth=min_rare_depth,
                          max_rare_depth=max_rare_depth,
                          status_update_callback=status_update_callback,
                          retain_intermediate_files=retain_intermediate_files)

if __name__ == "__main__":
    main()
