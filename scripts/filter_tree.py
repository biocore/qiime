#!/usr/bin/env python
# File created on 18 Jun 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = [
    "William Van Treuren",
    "Greg Caporaso",
    "Daniel McDonald",
    "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"


from cogent.parse.tree import DndParser
from qiime.filter import (filter_fasta, negate_tips_to_keep,
                          get_seqs_to_keep_lookup_from_seq_id_file,
                          get_seqs_to_keep_lookup_from_fasta_file,
                          filter_tree)
from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info[
    'brief_description'] = """This script prunes a tree based on a set of tip names"""

script_info[
    'script_description'] = """This script takes a tree and a list of OTU IDs (in one of several supported formats) and outputs a subtree retaining only the tips on the tree which are found in the inputted list of OTUs (or not found, if the --negate option is provided)."""

script_info['script_usage'] = []
script_info[
    'script_usage'].append(("""Prune a tree to include only the tips in tips_to_keep.txt""",
                            """""",
                            """%prog -i rep_seqs.tre -t tips_to_keep.txt -o pruned.tre"""))
script_info[
    'script_usage'].append(("""Prune a tree to remove the tips in tips_to_remove.txt. Note that the -n/--negate option must be passed for this functionality""",
                            """""",
                            """%prog -i rep_seqs.tre -t tips_to_keep.txt -o negated.tre -n"""))
script_info[
    'script_usage'].append(("""Prune a tree to include only the tips found in the fasta file provided""",
                            """""",
                            """%prog -i rep_seqs.tre -f fast_f.fna -o pruned_fast.tre"""))
script_info['output_description'] = \
    """Output is a pruned tree in newick format."""

script_info['required_options'] = [
    make_option('-i',
                '--input_tree_filepath',
                action='store',
                type='existing_filepath',
                dest='input_tree_fp',
                help='input tree filepath'),

    make_option('-o',
                '--output_tree_filepath',
                action='store',
                type='new_filepath',
                dest='output_tree_fp',
                help='output tree filepath'),
]

script_info['optional_options'] = [
    make_option('-n',
                '--negate',
                default=False,
                action='store_true',
                help='if negate is True will remove input tips/seqs, if \
   negate is False, will retain input tips/seqs [default: %default]'),

    make_option('-t',
                '--tips_fp',
                action='store',
                type='existing_filepath',
                help='A list of tips (one tip per line) or sequence identifiers \
  (tab-delimited lines with a seq identifier in the first field) \
  which should be retained \
  [default: %default]'),

    make_option('-f',
                '--fasta_fp',
                action='store',
                type='existing_filepath',
                help='A fasta file where the seq ids should be retained'
                ' [default: %default]'),
]

script_info['version'] = __version__


def main():

    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_tree_fp = opts.input_tree_fp
    tips_fp = opts.tips_fp
    fasta_fp = opts.fasta_fp
    output_tree_fp = opts.output_tree_fp

    if tips_fp is not None:
        tips_to_keep = get_seqs_to_keep_lookup_from_seq_id_file(
            open(tips_fp, 'U'))
    elif fasta_fp is not None:
        tips_to_keep = get_seqs_to_keep_lookup_from_fasta_file(
            open(fasta_fp, 'U'))
    else:
        option_parser.error("Must provide either -t or -f.")

    tree = DndParser(open(input_tree_fp, 'U'))

    if opts.negate:
        tips_to_keep = negate_tips_to_keep(tips_to_keep, tree)

    filtered_tree = filter_tree(tree, tips_to_keep)
    filtered_tree.writeToFile(output_tree_fp)

if __name__ == "__main__":
    # this comes in handy sometimes
    # import sys
    # sys.setrecursionlimit(10000)
    main()
