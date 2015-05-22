#!/usr/bin/env python
# File created on 10 Nov 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from qiime.util import parse_command_line_parameters, make_option
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode
from qiime.clean_raxml_parsimony_tree import decorate_numtips, decorate_depth,\
    get_insert_dict, drop_duplicate_nodes

scoring_methods = ['depth', 'numtips']

script_info = {}
script_info['brief_description'] = "Remove duplicate tips from Raxml Tree"
script_info[
    'script_description'] = "This script allows the user to remove specific duplicate tips from a Raxml tree."
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Example (depth):",
     "For this case the user can pass in input Raxml tree, duplicate tips, and define an output filepath. When using the depth option, only the deepest replicate is kept. ",
     " %prog -i raxml_v730_final_placement.tre -t 6 -o raxml_v730_final_placement_depth.tre"))
script_info['script_usage'].append(
    ("Example (numtips):",
     "For this case the user can pass in input Raxml tree, duplicate tips, and define an output filepath. When using the numtips option, the replicate with the fewest siblings is kept. ",
     " %prog -i raxml_v730_final_placement.tre -t 6 -o raxml_v730_final_placement_numtips.tre -s numtips"))
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option(
        '-i',
        '--input_tree',
        type="existing_filepath",
        help='the input raxml parsimony tree'),
    make_option(
        '-t',
        '--tips_to_keep',
        type="string",
        help='the input tips to score and retain (comma-separated list)'),
    make_option(
        '-o',
        '--output_fp',
        type="new_filepath",
        help='the output filepath'),
]
script_info['optional_options'] = [
    make_option(
        '-s',
        '--scoring_method',
        type="choice",
        help='the scoring method either depth or numtips [default: %default]',
        default='depth',
        choices=scoring_methods),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    # get options
    tree_fp = opts.input_tree
    tips_to_keep = opts.tips_to_keep.split(',')
    scoring_method = opts.scoring_method

    # load tree
    tree = DndParser(open(tree_fp, 'U'), constructor=PhyloNode)

    # decorate measurements onto tree (either by depth or by number of
    # children)
    if scoring_method == 'depth':
        tree2 = decorate_depth(tree)
    elif scoring_method == 'numtips':
        tree2 = decorate_numtips(tree)

    # get the nodes for the inserted sequences
    nodes_dict = get_insert_dict(tree2, set(tips_to_keep))

    # remove nodes accordingly
    final_tree = drop_duplicate_nodes(tree2, nodes_dict)

    # final_tree.nameUnnamedNodes()

    # write out the resulting tree
    open_outpath = open(opts.output_fp, 'w')
    open_outpath.write(final_tree.getNewick(with_distances=True))
    open_outpath.close()

if __name__ == "__main__":
    main()
