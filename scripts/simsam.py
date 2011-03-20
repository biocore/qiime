#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

import numpy
import random
import re 

from optparse import make_option

from cogent.parse.tree import DndParser

from qiime.parse import parse_otu_table
from qiime.util import parse_command_line_parameters, merge_n_otu_tables
from qiime.format import format_otu_table

from qiime.simsam import sort_nicely, get_new_otu_id, combine_sample_dicts

script_info = {}
script_info['brief_description'] = "Simulate samples for each sample in an OTU table, using a phylogenetic tree."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Make 3 related sample for each sample in in_otu_table.txt.","%prog -i in_otu_table.txt -t rep_set.tre -o out_otu_table.txt -d .001 -n 3")]
script_info['output_description']= "an otu table, samples are named: 'original_sample_0, original_sample_1 ...'"
script_info['required_options'] = [
 make_option('-i','--otu_table',help='the input otu table'),
 make_option('-t','--tree_file',help='tree file'),
 make_option('-o','--output_file',help='the output file'),
 make_option('-d','--dissim',help='dissimilarity between nodes up the tree',
    type=float),
 make_option('-n','--num',
    help='number of simulated samples per input sample',
    type=int)

]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
     parse_command_line_parameters(**script_info)

    out_fh = open(opts.output_file,'w')
    otu_table_fh = open(opts.otu_table,'U')
    sample_ids, otu_ids, otu_mtx, lineages = parse_otu_table(otu_table_fh)
    tree_fh = open(opts.tree_file,'U')
    tree = DndParser(tree_fh)
    tree_tips = [tip.Name for tip in tree.tips()]
    sort_nicely(tree_tips) # don't need sort here?
    # format of sample_dict: otu_id: num_seqs
    # sample_otu_tables = []

    # beware, get_new_otu_id may return something we already have in the table
    # don't want to rows of the same otu
    sample_dicts = []
    out_sam_names = []
    for i in range(len(sample_ids)):
        for j in range(opts.num):
            sample_dict = {}
            sample_vector = otu_mtx[:,i]
            for k in range(len(otu_ids)):
                if sample_vector[k] == 0: continue
                new_otu_id = get_new_otu_id(otu_ids[k], tree, opts.dissim)
                if sample_dict.has_key(new_otu_id):
                    sample_dict[new_otu_id] += sample_vector[k]
                else:
                    sample_dict[new_otu_id] = sample_vector[k]
            sample_dicts.append(sample_dict)
            out_sam_names.append(sample_ids[i] + '_' + str(j))


    res_otu_mtx, res_otus = combine_sample_dicts(sample_dicts)
    out_fh.write(format_otu_table(out_sam_names, res_otus, res_otu_mtx))


if __name__ == "__main__":
    main()