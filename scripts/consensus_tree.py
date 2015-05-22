#!/usr/bin/env python
# File created on 19 Aug 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "justinak@gmail.com"

import os

from skbio.tree import TreeNode, majority_rule

from qiime.parse import parse_newick, parse_rarefaction_fname
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                         make_option)


options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = "This script outputs a majority consensus tree given a collection of input trees."
script_info['script_description'] = ""
script_info['script_usage'] = [
    ("basic usage: given a directory of trees 'jackknifed_trees', compute the majority consensus and save as a newick formatted text file", "",
     "%prog -i jackknifed_trees -o consensus_tree.tre")

]

script_info[
    'output_description'] = "The output is a newick formatted tree compatible with most standard tree viewing programs"
script_info['required_options'] = [
    make_option('-i', '--input_dir', action='store',
                type='existing_dirpath', help='input folder containing trees'),
    make_option(
        '-o',
        '--output_fname',
        type='new_filepath',
        help='the output consensus tree filepath')
]
script_info['optional_options'] = [
    make_option('-s', '--strict',
                help='Use only nodes occurring >50% of the time [default: %default]',
                default=False, action='store_true')
]
script_info['version'] = __version__


def load_tree_files(tree_dir):
    """Load trees from filepaths

    checks if  filenames indicate that trees are from different
    distance methods.  If so, warns user.
    loads trees into phylonode objects
    returns [trees]
    raises a RuntimeError if no  trees are loaded
    """
    tree_file_names = os.listdir(tree_dir)
    # ignore invisible files like .DS_Store
    tree_file_names = [fname for fname in tree_file_names if not
                       fname.startswith('.')]

    # try to warn user if using multiple types of trees {
    try:
        base_names = []
        for fname in tree_file_names:
            base_names.append(parse_rarefaction_fname(fname)[0])
    except ValueError:
        pass
    else:
        if len(set(base_names)) > 1:
            warnstr = """
warning: trees are named differently, please be sure you're not
comparing trees generated in different manners, unless you're quite sure
that's what you intend to do.  types: """ + str(set(base_names)) + """
continuing anyway..."""
            warn(warnstr)
    # }
    trees = []
    for fname in tree_file_names:
        try:
            f = open(os.path.join(tree_dir, fname), 'U')
            tree = TreeNode.from_newick(f)
            tree.filepath = fname
            trees.append(tree)
            f.close()
        except IOError as err:
            sys.stderr.write('error loading tree ' + fname + '\n')
            exit(1)
    if len(trees) == 0:
        raise RuntimeError('Error: no trees loaded' +
                           ', check that tree directory has has valid trees')
    return trees


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    trees = load_tree_files(opts.input_dir)

    # this retains the default behavior from PyCogent's majority_rule function
    if opts.strict:
        cutoff = 0.5
    else:
        cutoff = 0.0
    consensus = majority_rule(trees=trees, cutoff=cutoff)[0]

    f = open(opts.output_fname, 'w')
    f.write(consensus.to_newick(with_distances=True))
    f.close()

if __name__ == "__main__":
    main()
