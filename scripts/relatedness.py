#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["William Van Treuren", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"
 
from cogent.util.option_parsing import (parse_command_line_parameters,
    make_option)
from biom.parse import parse_biom_table
from qiime.parse import parse_newick, PhyloNode
from qiime.relatedness_library import nri, nti

script_info = {}
script_info['brief_description'] = "Calculate NRI and NTI using formulas from Phylocom 4.2/3.41"
script_info['script_description'] = "This script calculates NRI and NTI from a path to a Newick formatted tree and a path to a comma separated list of ids in that tree that form the group whose NRI/NTI you want to test. The tree is not required to have distances. If none are found script will use the number of nodes (self inclusive) as their distance from one another. NRI and NTI are calculated as described in the Phylocom manual, not as in Webb 2002, or Webb 2000. The Phylocom manual is freely available on the web and Webb 2002 can be found in the Annual Review of Ecology and Systematics: Phylogenies and Community Ecology Webb 2002."
script_info['script_usage'] = [\
    ("Calculate both NRI and NTI from the given tree and group of taxa:",
     "",
     "%prog -t gg_tree.tre -i ids.txt -m nri,nti"),
    ("Calculate only NRI:",
     "",
     "%prog -t gg_tree.tre -i ids.txt -m nri"),
    ("Calculate only NTI using a different number of iterations:",
     "",
     "%prog -t gg_tree.tre -i ids.txt -m nti -i 100")]
script_info['output_description']= "Outputs a value for specified tests"
script_info['required_options'] = [\
 make_option('-t','--tree_fp',type="existing_filepath",help='the tree filepath'),
 make_option('-g','--taxa_fp',type="existing_filepath",help='taxa list filepath')]
script_info['optional_options'] = [\
 make_option('-i','--iters',type="int",default=1000,help='number of iterations to use for sampling tips without replacement (null model 2 community sampling, see see http://bodegaphylo.wikispot.org/Community_Phylogenetics). [default: %default]'),
 make_option('-m','--methods',type='string', default='nri,nti',help='comma-separated list of metrics to calculate. [default: %default]')]
script_info['version'] = __version__
script_info['help_on_no_arguments'] = True


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    tr = parse_newick(open(opts.tree_fp),PhyloNode)
    tip_dists, all_nodes = tr.tipToTipDistances() # tipTo returns a list of actual node objects
    all_ids = [node.Name for node in all_nodes]
    

    o = open(opts.taxa_fp)
    group_ids = [i.strip() for i in o.readline().split(',')]
    o.close()
    # check that there are at least 2 ids in the group, otherwise the math fails
    if len(group_ids) < 2:
        option_parser.error('you must have at least 2 taxa specified' +\
         ' in the taxa file or the math will fail.')

    # make sure specified taxa are in the tree, break at first failure
    for i in group_ids:
        try:
            all_ids.index(i)
        except ValueError:
            option_parser.error('taxa '+i+' not found in the tree')

    if len(all_ids)==len(group_ids): #m ust be the same set of ids if above check passes
        option_parser.error('The taxa_ids you specified contain every tip'+\
            ' in the tree. The NRI and NTI formulas will fail with these values'+\
            ' because there is no standard deviation of mpd or mntd, and thus'+\
            ' division by zero will occur. In addition, the concept of over/under'+\
            ' dispersion of a group of taxa (what NRI/NTI measure) is done in'+\
            ' reference to the tree they are a part of. If the group being tested'+\
            ' is the entire tree, the idea of over/under dispersion does not make'+\
            ' much sense.')

    # mapping from string of method name to function handle
    method_lookup = {'nri':nri, 'nti':nti}

    methods = opts.methods.split(',')
    for method in methods:
       if method not in method_lookup:
           option_parser.error("unknown method: %s; valid methods are: %s" % (method, ', '.join(method_lookup.keys())))
    
    for method in methods:
       print method+':', method_lookup[method](tip_dists, all_ids, group_ids, iters=opts.iters)

if __name__ == "__main__":
    main()
