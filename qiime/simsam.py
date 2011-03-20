#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

import numpy
import random
import re 

from optparse import make_option

from qiime.parse import parse_otu_table
from qiime.util import parse_command_line_parameters
from qiime.format import format_otu_table


#from codinghorror.com and Ned Batchelder
def sort_nicely( l ): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key ) 


def get_new_otu_id(old_otu_id, tree, dissim):
    node = tree.getNodeMatchingName(old_otu_id) #starts at tip
    uptree = 0
    while (not node.isRoot()) and uptree + node.Length < dissim:
        uptree += node.Length
        node = node.Parent

    # maybe do 50-50 descending each node, 
    # so we don't bias for overrepresented clades?
    if node.isTip():
        return node.Name
    else:
        return random.choice([tip.Name for tip in node.tips()])

def combine_sample_dicts(sample_dicts):
    all_otus = []
    for s in sample_dicts:
        all_otus.extend(s.keys())

    all_otus = list(set(all_otus))
    sort_nicely(all_otus)
    indices = {}
    for i in range(len(all_otus)):
        indices[all_otus[i]] = i

    otu_mtx = numpy.zeros((len(all_otus),len(sample_dicts)),int) # otu by sam
    for i, sample_dict in enumerate(sample_dicts):
        for otu, abund in sample_dict.items():
            otu_mtx[indices[otu],i] = abund

    return otu_mtx, all_otus