#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Jose Antonio Navas Molina",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from qiime.util import make_option, parse_command_line_parameters
from qiime.alpha_diversity import (single_file_alpha, multiple_file_alpha,
                                   list_known_metrics)
import os

script_info = {}
script_info['brief_description'] = """Calculate alpha diversity on each sample in an otu table, using a variety of alpha diversity metrics"""
script_info['script_description'] = \
"""This script calculates alpha diversity, or within-sample diversity, using an
OTU table. The QIIME pipeline allows users to conveniently calculate more than
two dozen different diversity metrics. The full list of available metrics is
available by passing the -s option to this script.

Documentation of the metrics can be found at
http://scikit-bio.org/docs/latest/generated/skbio.diversity.alpha.html.
Every metric has different strengths and limitations - technical discussion of
each metric is readily available online and in ecology textbooks, but is beyond
the scope of this document.
"""
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Single File Alpha Diversity Example (non-phylogenetic):""",
     """To perform alpha diversity (e.g. chao1) on a single OTU table, where the results are output to "alpha_div.txt", you can use the following command:""",
     """%prog -i otu_table.biom -m chao1 -o adiv_chao1.txt"""))

script_info['script_usage'].append(
    ("""Single File Alpha Diversity Example (phylogenetic):""",
     """In the case that you would like to perform alpha diversity using a phylogenetic metric (e.g. PD_whole_tree), you can use the following command:""",
     """%prog -i otu_table.biom -m PD_whole_tree -o adiv_pd.txt -t rep_set.tre"""))

script_info['script_usage'].append(
    ("""Single File Alpha Diversity Example with multiple metrics:""",
     """You can use the following idiom to run multiple metrics at once (comma-separated):""",
     """%prog -i otu_table.biom -m chao1,PD_whole_tree -o adiv_chao1_pd.txt -t rep_set.tre"""))

script_info['script_usage'].append(
    ("""Multiple File (batch) Alpha Diversity:""",
     """To perform alpha diversity on multiple OTU tables (e.g.: rarefied otu tables resulting from multiple_rarefactions.py), specify an input directory instead of a single otu table, and an output directory (e.g. "alpha_div_chao1_PD/") as shown by the following command:""",
     """%prog -i otu_tables/ -m chao1,PD_whole_tree -o adiv_chao1_pd/ -t rep_set.tre"""))

script_info['output_description'] = """The resulting file(s) is a tab-delimited text file, where the columns correspond to alpha diversity metrics and the rows correspond to samples and their calculated diversity measurements. When a folder is given as input (-i), the script processes every otu table file in the given folder, and creates a corresponding file in the output directory.

Example Output:

====== ======= ============= =============
\      simpson PD_whole_tree observed_otus
====== ======= ============= =============
PC.354 0.925   2.83739       16.0
PC.355 0.915   3.06609       14.0
PC.356 0.945   3.10489       19.0
PC.481 0.945   3.65695       19.0
PC.593 0.91    3.3776        15.0
PC.607 0.92    4.13397       16.0
PC.634 0.9     3.71369       14.0
PC.635 0.94    4.20239       18.0
PC.636 0.925   3.78882       16.0
====== ======= ============= =============
"""
script_info['required_options'] = []
script_info['optional_options'] = [
    make_option('-i', '--input_path',
                help='Input OTU table filepath or input directory containing OTU' +
                ' tables for batch processing. [default: %default]',
                type='existing_path'),
    make_option('-o', '--output_path',
                help='Output filepath to store alpha diversity metric(s) for each sample in a tab-separated format '
                'or output directory when batch processing. [default: %default]',
                type='new_path'),
    make_option('-m', '--metrics', type='multiple_choice',
                mchoices=list_known_metrics(),
                default='PD_whole_tree,chao1,observed_otus',
                help='Alpha-diversity metric(s) to use. A comma-separated list should' +
                ' be provided when multiple metrics are specified. [default: %default]'),
    make_option('-s', '--show_metrics', action='store_true',
                dest="show_metrics",
                help='Show the available alpha-diversity metrics and exit.'),
    make_option('-t', '--tree_path', default=None,
                help='Input newick tree filepath.' +
                ' [default: %default; REQUIRED for phylogenetic metrics]',
                type='existing_filepath')
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    if opts.show_metrics:
        print("Known metrics are: %s\n"
              % (', '.join(list_known_metrics()),))
        print("For more information, see http://scikit-bio.org/docs/latest/"
              "generated/skbio.diversity.alpha.html")
        exit(0)
    almost_required_options = ['input_path', 'output_path', 'metrics']
    for option in almost_required_options:
        if getattr(opts, option) is None:
            option_parser.error('Required option --%s omitted.' % option)

    if os.path.isdir(opts.input_path):
        multiple_file_alpha(opts.input_path, opts.output_path, opts.metrics,
                            opts.tree_path)
    elif os.path.isfile(opts.input_path):
        try:
            f = open(opts.output_path, 'w')
            f.close()
        except IOError:
            if os.path.isdir(opts.output_path):
                option_parser.error(
                    "ioerror, couldn't create output file. The output path is a directory, which should be a single file")
            else:
                option_parser.error("ioerror, couldn't create output file")
        single_file_alpha(opts.input_path, opts.metrics,
                          opts.output_path, opts.tree_path)


if __name__ == "__main__":
    main()
