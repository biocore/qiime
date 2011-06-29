#!/usr/bin/env python
# File created on 18 Jun 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren","Greg Caporaso", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"
__status__ = "Release"


from cogent.parse.tree import DndParser
from qiime.filter import (filter_fasta, negate_tips_to_keep,
                          get_seqs_to_keep_lookup_from_seq_id_file,
                          get_seqs_to_keep_lookup_from_fasta_file)
from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = """This script prunes a tree based on a set of tip names"""
                                      
script_info['script_description'] = """This script takes a tree and a list of OTU IDs (in one of several supported formats) and outputs a subtree retaining only the tips on the tree which are found in the inputted list of OTUs (or not found, if the --negate option is provided)."""
    
script_info['script_usage'] = []
script_info['script_usage'].append(("""Prune a tree to include only the tips in tips_to_keep.txt""",\
    """""",\
    """%prog -i rep_seqs.tre -t tips_to_keep.txt -o rep_seqs_subtree.tre"""))
script_info['script_usage'].append(("""Prune a tree to remove the tips in tips_to_remove.txt. Note that the -n/--negate option must be passed for this functionality.""",\
    """""",\
    """%prog -i rep_seqs.tre -t tips_to_remove.txt -o rep_seqs_subtree.tre -n"""))

script_info['output_description'] = \
    """Output is a pruned tree in newick format."""

script_info['required_options']=[\
    make_option('-i','--input_tree_fp',help='input tree filepath'),
    make_option('-o','--output_tree_fp',help='output tree filepath')
]

script_info['optional_options']=[\
    make_option('-n','--negate',default=False, action='store_true',
            help='if negate is not false will prune tips fed in and save \
            all others [default: %default]'),
    make_option('-t','--tips_fp', 
            help='A list of sequence identifiers (or tab-delimited lines with'
                 ' a seq identifier in the first field) which should be retained'
                 ' [default: %default]'),
    make_option('-f','--fasta_fp',
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
    
    if tips_fp != None:
        tips_to_keep = get_seqs_to_keep_lookup_from_seq_id_file(open(tips_fp,'U'))
    elif fasta_fp != None:
        tips_to_keep = get_seqs_to_keep_lookup_from_fasta_file(open(fasta_fp,'U'))
    else:
        option_parser.error("Must provide either -t or -f.")
    
    tree = DndParser(open(input_tree_fp,'U'))
    
    if opts.negate:
        tips_to_keep = negate_tips_to_keep(tips_to_keep, tree)
    
    tree_out = tree.getSubTree(tips_to_keep)
   
    tree_out.writeToFile(output_tree_fp)

if __name__ == "__main__":
    main()
