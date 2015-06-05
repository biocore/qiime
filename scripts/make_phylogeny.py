#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski, Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from os.path import splitext
from qiime.util import make_option

from qiime.util import parse_command_line_parameters

from qiime.make_phylogeny import tree_module_names, tree_method_constructors,\
    CogentTreeBuilder

import warnings

script_info = {}
script_info['brief_description'] = """Make Phylogeny"""
script_info[
    'script_description'] = """Many downstream analyses require that the phylogenetic tree relating the OTUs in a study be present. The script make_phylogeny.py produces this tree from a multiple sequence alignment. Trees are constructed with a set of sequences representative of the OTUs, by default using FastTree (Price, Dehal, & Arkin, 2009)."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Examples:""",
     """A simple example of make_phylogeny.py is shown by the following command, where we use the default tree building method (fasttree) and write the file to the current working directory without a log file:""",
     """%prog -i $PWD/aligned.fasta -o $PWD/rep_phylo.tre"""))
script_info['script_usage'].append(
    ("""""",
     """Alternatively, if the user would prefer using another tree building method (i.e. clearcut (Sheneman, Evans, & Foster, 2006)), then they could use the following command:""",
     """%prog -i $PWD/aligned.fasta -t clearcut"""))
script_info['output_description'] = """The result of make_phylogeny.py consists of a newick formatted tree file (.tre) and optionally a log file. The tree file is formatted using the Newick format and this file can be viewed using most tree visualization tools, such as TopiaryTool, FigTree, etc.

The tips of the tree are the first word from the input sequences from the fasta file, e.g.: '>101 PC.481_71 RC:1..220' is represented in the tree as '101'."""
script_info['required_options'] = [
    make_option('-i', '--input_fp', action='store',
                type='existing_filepath', dest='input_fp', help='Path to read ' +
                'input fasta alignment, only first word in defline will be considered')
]
valid_root_methods = ['midpoint', 'tree_method_default']

script_info['optional_options'] = [
    make_option(
        '-t', '--tree_method', action='store', type='choice', choices=list(tree_module_names.keys()),
        help='Method for tree building. Valid choices are: ' +
        ', '.join(tree_module_names.keys()) +
        ' [default: %default]', default='fasttree'),

    make_option('-o', '--result_fp', action='store', type='new_filepath',
                help='Path to store ' +
                'result file [default: <input_sequences_filename>.tre]'),

    make_option('-l', '--log_fp', action='store', type='new_filepath',
                help='Path to store ' +
                'log file [default: No log file created.]'),

    make_option(
        '-r', '--root_method', action='store', type='choice', choices=list(valid_root_methods),
        help='method for choosing root of phylo tree' +
        '  Valid choices are: ' + ', '.join(valid_root_methods) +
        ' [default: tree_method_default]',
        default='tree_method_default'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if not (opts.tree_method in tree_method_constructors or
            opts.tree_method in tree_module_names):
        option_parser.error(
            'Invalid alignment method: %s.\nValid choices are: %s'
            % (opts.tree_method,
               ' '.join(tree_method_constructors.keys() +
                        tree_module_names.keys())))
    try:
        tree_builder_constructor =\
            tree_method_constructors[opts.tree_method]
        tree_builder_type = 'Constructor'
        params = {}
        tree_builder = tree_builder_constructor(params)
    except KeyError:
        tree_builder = CogentTreeBuilder({
            'Module': tree_module_names[opts.tree_method],
            'Method': opts.tree_method
        })
        tree_builder_type = 'Cogent'

    input_seqs_filepath = opts.input_fp
    result_path = opts.result_fp
    if not result_path:  # empty or None
        fpath, ext = splitext(input_seqs_filepath)  # fpath omits extension
        result_path = fpath + ".tre"

    open(result_path, 'w').close()  # touch
    log_path = opts.log_fp
    if log_path is not None:
        open(log_path, 'w').close()

    if tree_builder_type == 'Constructor':
        tree_builder(input_seqs_filepath,
                     result_path=result_path, log_path=log_path, failure_path=failure_path)
    elif tree_builder_type == 'Cogent':
        tree_builder(result_path, aln_path=input_seqs_filepath,
                     log_path=log_path, root_method=opts.root_method)


if __name__ == "__main__":
    main()
