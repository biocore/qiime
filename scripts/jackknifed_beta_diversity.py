#!/usr/bin/env python
# File created on 19 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = [
    "Greg Caporaso",
    "Justin Kuczynski",
    "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import make_option
from os import makedirs
from qiime.util import (load_qiime_config,
                        parse_command_line_parameters,
                        get_options_lookup)
from qiime.parse import parse_qiime_parameters
from qiime.workflow.util import (print_commands,
                                 call_commands_serially,
                                 print_to_stdout,
                                 no_status_updates,
                                 validate_and_set_jobs_to_start)
from qiime.workflow.downstream import run_jackknifed_beta_diversity

script_info = {}

script_info['brief_description'] = """A workflow script for performing jackknifed\
 UPGMA clustering and building jackknifed Emperor PCoA plots."""

script_info['script_description'] = """To directly measure the robustness of\
 individual UPGMA clusters and clusters in PCoA plots, one can\
 perform jackknifing (repeatedly resampling a subset of the available data\
 from each sample)."""

script_info['script_usage'] = []

script_info['script_usage'].append(("""Example:""", """These steps are performed\
 by the following command: Compute beta diversity distance matrix from otu\
 table (and tree, if applicable); build rarefied OTU tables by evenly sampling\
 to the specified depth (-e); build UPGMA tree from full distance matrix;\
 compute distance matrics for rarefied OTU tables; build UPGMA trees from\
 rarefied OTU table distance matrices; build a consensus tree from the rarefied\
 UPGMA trees; compare rarefied OTU table distance matrix UPGMA trees to either\
 (full or consensus) tree for jackknife support of tree nodes; perform\
 principal coordinates analysis on distance matrices generated from rarefied\
 OTU tables; generate Emperor PCoA plots with jackknifed support.

""", """%prog -i otu_table.biom -o bdiv_jk100 -e 100 -m Fasting_Map.txt\
 -t rep_set.tre"""))

script_info['script_usage_output_to_remove'] = ['bdiv_jk100']

script_info['output_description'] = """This scripts results in several distance\
 matrices (from beta_diversity.py), several rarefied OTU tables\
 (from multiple_rarefactions_even_depth.py), several UPGMA trees (from upgma_cluster.py),\
 a supporting file and newick tree with support values (from tree_compare.py),\
 and Emperor PCoA plots."""

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info['required_options'] = [
    make_option('-i', '--otu_table_fp', type='existing_filepath',
                help='the input OTU table in biom format [REQUIRED]'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output directory [REQUIRED]'),
    make_option('-e', '--seqs_per_sample', type='int',
                help='number of sequences to include in each jackknifed subset' +
                ' [REQUIRED]'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='path to the mapping file [REQUIRED]'),
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
    make_option('--master_tree', default="consensus",
                type='choice', choices=['consensus', 'full'],
                help='method for computing master trees in jackknife analysis.' +
                ' "consensus": consensus of trees from jackknifed otu tables. ' +
                ' "full": tree generated from input (unsubsambled) otu table. ' +
                ' [default: %default]'),
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
    options_lookup['jobs_to_start_workflow']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = \
        parse_command_line_parameters(**script_info)

    verbose = opts.verbose

    otu_table_fp = opts.otu_table_fp
    output_dir = opts.output_dir
    tree_fp = opts.tree_fp
    seqs_per_sample = opts.seqs_per_sample
    verbose = opts.verbose
    print_only = opts.print_only
    master_tree = opts.master_tree

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

    run_jackknifed_beta_diversity(otu_table_fp=otu_table_fp,
                                  tree_fp=tree_fp,
                                  seqs_per_sample=seqs_per_sample,
                                  output_dir=output_dir,
                                  command_handler=command_handler,
                                  params=params,
                                  qiime_config=qiime_config,
                                  mapping_fp=opts.mapping_fp,
                                  parallel=parallel,
                                  status_update_callback=status_update_callback,
                                  master_tree=master_tree)


if __name__ == "__main__":
    main()
