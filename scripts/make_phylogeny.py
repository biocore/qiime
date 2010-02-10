#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski, Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 
from os.path import splitext
from optparse import make_option

from qiime.util import parse_command_line_parameters

from qiime.make_phylogeny import tree_module_names, tree_method_constructors,\
    CogentTreeBuilder

script_description = """Build a phylogenetic tree from aligned sequences, using one of several techniques."""

script_usage = """
 Build phylogenetic tree from sequences in  align/rep_set_aligned.fa
 writing tree to align/rep_set_aligned.tre using FastTree (default):
    make_phylogeny.py -i align/rep_set_aligned.fa -o align/rep_set_aligned.tre -l align/rep_set_tree.log
"""

required_options = [
    make_option('-i','--input_fp',action='store',
          type='string',dest='input_fp',help='Path to read '+\
          'input alignment')
]

optional_options = [\
    make_option('-t','--tree_method',action='store',
          help='Method for tree building. Valid choices are: '+\
          ', '.join(tree_module_names.keys())+\
          ' [default: %default]', default='fasttree'),
          
    make_option('-o','--result_fp',action='store',
          help='Path to store '+\
          'result file [default: <input_sequences_filename>.tre]'),
          
    make_option('-l','--log_fp',action='store',
          help='Path to store '+\
          'log file [default: No log file created.]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if not (opts.tree_method in tree_method_constructors or
            opts.tree_method in tree_module_names):
        option_parser.error(\
         'Invalid alignment method: %s.\nValid choices are: %s'\
         % (opts.tree_method,\
            ' '.join(tree_method_constructors.keys() +
                tree_module_names.keys())))
    
    #verbose = opts.verbose
    try:
        tree_builder_constructor =\
            tree_method_constructors[opts.tree_method]
        tree_builder_type = 'Constructor'
        params = {}
        tree_builder = tree_builder_constructor(params)
    except KeyError:
        tree_builder = CogentTreeBuilder({
                'Module':tree_module_names[opts.tree_method],
                'Method':opts.tree_method
                })
        tree_builder_type = 'Cogent'
     
    input_seqs_filepath = opts.input_fp
    result_path = opts.result_fp
    if not result_path: # empty or None
        fpath, ext = splitext(input_seqs_filepath) # fpath omits extension
        result_path = fpath + ".tre"
     
    log_path = opts.log_fp

    if tree_builder_type=='Constructor':
        tree_builder(input_seqs_filepath,\
        result_path=result_path,log_path=log_path,failure_path=failure_path)
    elif tree_builder_type=='Cogent':
        tree_builder(result_path, aln_path=input_seqs_filepath,
            log_path=log_path)


if __name__ == "__main__":
    main()
