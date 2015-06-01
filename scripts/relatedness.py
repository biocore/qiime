#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["William Van Treuren", "Jose Carlos Clemente Litran",
               "Jose Antonio Navas Molina", "Yoshiki Vazquez Baeza",
               "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"

from sys import stdout

from qiime.parse import parse_newick, PhyloNode
from qiime.relatedness_library import nri, nti
from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = (
    "Calculate NRI (net relatedness index) and NTI (nearest taxon index) "
    "using the formulas from Phylocom 4.2/3.41 and Webb 2002.")
script_info['script_description'] = (
    "This script calculates NRI and NTI from a path to a Newick formatted "
    "tree and a path to a comma separated list of ids in that tree that form "
    "the group whose NRI/NTI you want to test. The tree is not required to "
    "have distances. If none are found script will use the number of nodes "
    "(self inclusive) as their distance from one another. NRI and NTI are "
    "calculated as described in the Phylocom manual (which is a slightly "
    "modified version of that found in Webb 2002, and Webb 2000). The "
    "Phylocom manual is freely available on the web and Webb 2002 can be "
    "found in the Annual Review of Ecology and Systematics: Phylogenies and "
    "Community Ecology Webb 2002.")
script_info['script_usage'] = [
    ("Calculate both NRI and NTI from the given tree and group of taxa:",
     "",
     "%prog -t reference.tre -g group1_otus.txt -m nri,nti"),
    ("Calculate only NRI:",
     "",
     "%prog -t reference.tre -g group1_otus.txt -m nri"),
    ("Calculate only NTI using a different number of iterations:",
     "",
     "%prog -t reference.tre -g group1_otus.txt -m nti -i 100"),
    ("Calculate only NTI using a different number of iterations and save the "
     "results into a file called output.txt", "", "%prog -t reference.tre -g "
     "group1_otus.txt -m nti -i 100 -o output.txt")]
script_info['output_description'] = "Outputs a value for specified tests"
script_info['required_options'] = [
    make_option('-t', '--tree_fp', type="existing_filepath",
                help='the tree filepath'),
    make_option('-g', '--taxa_fp', type="existing_filepath",
                help='taxa list filepath')]
script_info['optional_options'] = [
    make_option('-i', '--iters', type="int", default=1000,
                help='number of iterations to use for sampling tips without '
                     'replacement (null model 2 community sampling, see '
                     'http://bodegaphylo.wikispot.org/Community_Phylogenetics '
                     '[default: %default]'),
    make_option('-m', '--methods', type='multiple_choice',
                default='nri,nti', mchoices=['nri', 'nti'],
                help='comma-separated list of metrics to calculate. '
                     '[default: %default]'),
    make_option('-o', '--output_fp', type="new_filepath",
                help="path where output will be written [default: print to "
                     "screen]", default=None)]
script_info['version'] = __version__
script_info['help_on_no_arguments'] = True


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_fp = opts.output_fp
    if output_fp:
        fd = open(output_fp, 'w')
    else:
        fd = stdout

    tr = parse_newick(open(opts.tree_fp), PhyloNode)
    # all_nodes is list node objs
    tip_dists, all_nodes = tr.tipToTipDistances()
    all_ids = [node.Name for node in all_nodes]

    o = open(opts.taxa_fp)
    group_ids = [i.strip() for i in o.readline().split(',')]
    o.close()
    # check that there are at least 2 ids in the group, otherwise the math
    # fails
    if len(group_ids) < 2:
        option_parser.error('Not enough taxa in the taxa file.You must have '
                            'at least 2 taxa specified in the taxa file or '
                            'the standard deviation of the distance will be '
                            'zero, causing both NRI and NTI to fail.')
    # check that all_ids contains every group_id
    if not set(group_ids).issubset(all_ids):
        raise option_parser.error('There are taxa in the taxa file which are '
                                  'not found in the tree. You may have '
                                  'specified an internal node.')
    # check that all_ids != group_ids
    if len(all_ids) == len(group_ids):  # must be same set if above passes
        option_parser.error('The taxa_ids you specified contain every tip in '
                            'the tree. The NRI and NTI formulas will fail '
                            'because there is no standard deviation of mpd or '
                            'mntd, and thus division by zero. In addition, '
                            'the concept of over/under dispersion of a group '
                            'of taxa (what NRI/NTI measure) is done in '
                            'reference to the tree they are a part of. If the '
                            'group being tested is the entire tree, the idea '
                            'of over/under dispersion makes little sense.')

    # mapping from string of method name to function handle
    method_lookup = {'nri': nri, 'nti': nti}

    methods = opts.methods
    for method in methods:
        if method not in method_lookup:
            option_parser.error("Unknown method: %s; valid methods are: %s" %
                                (method, ', '.join(method_lookup.keys())))

    for method in methods:
        print >> fd, method + ':', method_lookup[method](tip_dists, all_ids,
                                                         group_ids,
                                                         iters=opts.iters)

    fd.close()

if __name__ == "__main__":
    main()
