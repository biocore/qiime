#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

from cogent.parse.tree import DndParser
from qiime.util import parse_command_line_parameters, make_option
from qiime.simsam import sim_otu_table
from qiime.format import format_biom_table
from biom.table import table_factory
from biom.parse import parse_biom_table

script_info = {}
script_info['brief_description'] = "Simulate samples for each sample in an OTU table, using a phylogenetic tree."
script_info['script_description'] = """ This script makes n samples related to each sample in an input otu table

An input OTU table with 3 samples and n=2 will result in an output otu table with 6 samples total: 3 clusters of 2 related samples.

To simulate each of new samples, this script uses a sample in the input OTU table, and for each OTU in that sample the script
traverses rootward on the tree a distance specified by '-d' to a point x. It then randomly selects a tip that decends from x,
(call that new tip 'o2'), and reassigns all observations of the original OTU to the tip/OTU 'o2'.
"""
script_info['script_usage'] = [("","Make 3 related sample for each sample in otu_table.biom.","%prog -i otu_table.biom -t rep_set.tre -o otu_table.simsam.biom -d .001 -n 3")]
script_info['output_description']= "an otu table, samples are named: 'original_sample_0, original_sample_1 ...'"
script_info['required_options'] = [
 make_option('-i','--otu_table',help='the input otu table',type='existing_filepath'),
 make_option('-t','--tree_file',help='tree file',type='existing_filepath'),
 make_option('-o','--output_file',help='the output file',type='new_filepath'),
 make_option('-d','--dissim',help='dissimilarity between nodes up the tree',
    type='float'),
 make_option('-n','--num',
    help='number of simulated samples per input sample',
    type='int')

]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
     parse_command_line_parameters(**script_info)

    out_fh = open(opts.output_file,'w')
    otu_table_fh = open(opts.otu_table,'U')
    otu_table = parse_biom_table(otu_table_fh)
    tree_fh = open(opts.tree_file,'U')
    tree = DndParser(tree_fh)

    res_sam_names, res_otus, res_otu_mtx, res_otu_metadata = \
     sim_otu_table(otu_table.SampleIds, otu_table.ObservationIds, otu_table.iterSamples(), 
                   otu_table.ObservationMetadata, tree, opts.num, opts.dissim)


    rich_table = table_factory(res_otu_mtx,res_sam_names,res_otus,
    observation_metadata=res_otu_metadata)
    out_fh.write(format_biom_table(rich_table))


if __name__ == "__main__":
    main()
