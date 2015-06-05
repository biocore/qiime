#!/usr/bin/env python
# File created on 04 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Kestrel Gorlick"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from skbio.util import create_dir

from qiime.util import make_option
from os import makedirs
from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup
from qiime.parse import parse_qiime_parameters
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 print_to_stdout,
                                 no_status_updates,
                                 validate_and_set_jobs_to_start)
from qiime.workflow.downstream import run_beta_diversity_through_plots


# beta_diversity_through_3d_plots.py
qiime_config = load_qiime_config()
options_lookup = get_options_lookup()
script_info = {}
script_info[
    'brief_description'] = """A workflow script for computing beta diversity distance matrices and generating PCoA plots"""
script_info['script_description'] = """This script will perform beta diversity, principal coordinate analysis, and generate a preferences file along with 3D PCoA Plots.
"""
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Example:""",
     """Given an OTU table, a phylogenetic tree, an even sampling depth, and a mapping file, perform the following steps: 1. Randomly subsample otu_table.biom to even number of sequences per sample (100 in this case); 2. Compute a weighted and unweighted unifrac distance matrices (can add additional metrics by passing a parameters file via -p); 3. Peform a principal coordinates analysis on the result of Step 2; 4. Generate a 2D and 3D plots for all mapping fields.""",
     """%prog -i otu_table.biom -o bdiv_even100/ -t rep_set.tre -m Fasting_Map.txt -e 100"""))

script_info['script_usage_output_to_remove'] = ['bdiv_even100']

script_info[
    'output_description'] = """This script results in a distance matrix (from beta_diversity.py), a principal coordinates file (from principal_coordinates.py), and folders containing the resulting PCoA plots (accessible through html files)."""

script_info['required_options'] = [
    make_option('-i', '--otu_table_fp', type='existing_filepath',
                help='the input biom table [REQUIRED]'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='path to the mapping file [REQUIRED]'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output directory [REQUIRED]'),
]
script_info['optional_options'] = [
    make_option('-t', '--tree_fp', type='existing_filepath',
                help='path to the tree file [default: %default; ' +
                'REQUIRED for phylogenetic measures]'),
    make_option('-p', '--parameter_fp', type='existing_filepath',
                help='path to the parameter file, which specifies changes' +
                ' to the default behavior. ' +
                'See http://www.qiime.org/documentation/file_formats.html#qiime-parameters .' +
                ' [if omitted, default values will be used]'),
    make_option('--color_by_all_fields',
                help='plots will have coloring for all mapping fields ' +
                '[default: %default; only include fields with greater than one value ' +
                'and fewer values than the number of samples]',
                default=False, action='store_true'),
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
    make_option('-e', '--seqs_per_sample', type='int',
                help='depth of coverage for even sampling [default: %default]'),
    make_option('--suppress_emperor_plots', action='store_true',
                help='Do not generate emperor plots [default: %default]',
                default=False),
    options_lookup['jobs_to_start_workflow']]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose

    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    verbose = opts.verbose
    print_only = opts.print_only
    seqs_per_sample = opts.seqs_per_sample

    parallel = opts.parallel
    # No longer checking that jobs_to_start > 2, but
    # commenting as we may change our minds about this.
    #if parallel: raise_error_on_parallel_unavailable()

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

    create_dir(output_dir, fail_on_exist=not opts.force)

    if print_only:
        command_handler = print_commands
    else:
        command_handler = call_commands_serially

    if verbose:
        status_update_callback = print_to_stdout
    else:
        status_update_callback = no_status_updates

    run_beta_diversity_through_plots(otu_table_fp=otu_table_fp,
                                     mapping_fp=mapping_fp,
                                     output_dir=output_dir,
                                     command_handler=command_handler,
                                     params=params,
                                     qiime_config=qiime_config,
                                     color_by_interesting_fields_only=not opts.color_by_all_fields,
                                     sampling_depth=seqs_per_sample,
                                     tree_fp=tree_fp,
                                     parallel=parallel,
                                     suppress_emperor_plots=opts.suppress_emperor_plots,
                                     status_update_callback=status_update_callback)

if __name__ == "__main__":
    main()
