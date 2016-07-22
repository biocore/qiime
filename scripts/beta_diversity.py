#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from qiime.util import make_option, parse_command_line_parameters
from qiime.beta_diversity import (single_file_beta, multiple_file_beta,
                                  list_known_metrics)
import os
from sys import stderr

script_info = {}
script_info['brief_description'] = """Calculate beta diversity (pairwise sample\
 dissimilarity) on one or many otu tables"""
script_info['script_description'] = """The input for this script is the OTU table\
 containing the number of sequences observed in each OTU (rows) for each sample\
 (columns). For more information pertaining to the OTU table refer to the\
 documentation for make_otu_table. If the user would like phylogenetic beta\
 diversity metrics using UniFrac, a phylogenetic tree must also be passed as\
 input (see make_phylogeny.py). The output of this script is a distance matrix\
 containing a dissimilarity value for each pairwise comparison.

A number of metrics are currently supported, including unweighted and weighted\
 UniFrac (pass the -s option to see available metrics). In general, because\
 unifrac uses phylogenetic information, one of the unifrac metrics is\
 recommended, as results can be vastly more useful (Hamady & Knight, 2009).\
 Quantitative measures (e.g. weighted unifrac) are ideally suited to revealing\
 community differences that are due to changes in relative taxon abundance\
 (e.g., when a particular set of taxa flourish because a limiting nutrient\
 source becomes abundant). Qualitative measures (e.g. unweighted unifrac) are\
 most informative when communities differ primarily by what can live in them\
 (e.g., at high temperatures), in part because abundance information can\
 obscure significant patterns of variation in which taxa are present (Lozupone\
 et al., 2007). Most qualitative measures are referred to here e.g.\
 "binary_jaccard". Typically both weighted and unweighted unifrac are used."""
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Single File Beta Diversity (non-phylogenetic):""",
     """To perform beta diversity (using e.g. euclidean distance) on a single OTU\
 table, where the results are output to beta_div/, use the following command:""",
     """%prog -i otu_table.biom -m euclidean -o beta_div"""))

script_info['script_usage'].append(
    ("""Single File Beta Diversity (phylogenetic):""",
     """In the case that you would like to perform beta diversity using a\
 phylogenetic metric (e.g. weighted_unifrac), you can use the following\
 command:""",
     """%prog -i otu_table.biom -m weighted_unifrac -o beta_div/ -t rep_set.tre"""))

script_info['script_usage'].append(
    ("""Multiple File (batch) Beta Diversity (phylogenetic):""",
     """To perform beta diversity on multiple OTU tables (e.g., resulting files from\
 multiple_rarefactions.py), specify an input directory (e.g. otu_tables/) as\
 shown by the following command:""",
     """%prog -i otu_tables/ -m weighted_unifrac -o beta_div/ -t rep_set.tre"""))

script_info['output_description'] = """Each file in the input directory should be\
 an otu table, and the output of beta_diversity.py is a folder containing text\
 files, each a distance matrix between samples corresponding to an input otu\
 table."""
script_info['required_options'] = []
script_info['optional_options'] = [
    make_option('-i', '--input_path',
                help='Input OTU table in biom format or input directory containing OTU ' +
                'tables in biom format for batch processing.',
                type='existing_path'),
    make_option('-r', '--rows', default=None, type='string',
                help='Compute for only these rows of the distance matrix.' +
                ' User should pass a list of sample names (e.g. "s1,s3")' +
                ' [default: %default; full n x n matrix is generated]'),
    make_option('-o', '--output_dir',
                help="Output directory. One will be created if it doesn't exist.",
                type='new_dirpath'),
    make_option(
        '-m', '--metrics', default='unweighted_unifrac,weighted_unifrac',
        type='multiple_choice', mchoices=list_known_metrics(),
        help='Beta-diversity metric(s) to use. A comma-separated list should' +
        ' be provided when multiple metrics are specified. [default: %default]'),
    make_option('-s', '--show_metrics', action='store_true', default=False,
                help='Show the available beta-diversity metrics and exit. Metrics' +
                ' starting with' +
                ' "binary..." specifies that a metric is qualitative, and considers' +
                ' only the presence or absence of each taxon [default: %default]'),
    make_option('-t', '--tree_path', default=None,
                help='Input newick tree filepath, which is required when phylogenetic' +
                ' metrics are specified. [default: %default]',
                type='existing_filepath'),
    make_option('-f', '--full_tree', action="store_true", default=False,
                help='By default, tips not corresponding to OTUs in the OTU table are ' +
                'removed from the tree for diversity calculations. ' +
                'Pass to skip this step if you\'re already passing a minimal tree.' +
                ' Beware with "full_tree" metrics, as extra tips in the tree' +
                ' change the result'),
]
script_info['option_label'] = {'input_path': 'OTU table filepath',
                               'rows': 'List of samples for compute',
                               'metrics': 'Metrics to use',
                               'show_metrics': 'Show metrics',
                               'tree_path': 'Newick tree filepath',
                               'full_tree': 'Tree already trimmed',
                               'output_dir': 'Output directory'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.show_metrics:
        print("Known metrics are: %s\n"
              % (', '.join(list_known_metrics()),))
        exit(0)

    almost_required_options = ['input_path', 'output_dir', 'metrics']
    for option in almost_required_options:
        if getattr(opts, option) is None:
            option_parser.error('Required option --%s omitted.' % option)

    if opts.output_dir.endswith('.txt'):
        stderr.write('output must be a directory, files will be named' +
                     ' automatically.  And we refuse to make .txt directories\n')
        exit(1)

    if opts.tree_path == "None":
        opts.tree_path = None

    try:
        os.makedirs(opts.output_dir)
    except OSError:
        pass  # hopefully dir already exists

    if os.path.isdir(opts.input_path):
        multiple_file_beta(opts.input_path, opts.output_dir, opts.metrics,
                           opts.tree_path, opts.rows, full_tree=opts.full_tree)
    elif os.path.isfile(opts.input_path):
        single_file_beta(opts.input_path, opts.metrics, opts.tree_path,
                         opts.output_dir, opts.rows, full_tree=opts.full_tree)
    else:
        stderr.write("io error, input path not valid.  Does it exist?")
        exit(1)


if __name__ == "__main__":
    main()
