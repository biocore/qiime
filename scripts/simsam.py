#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

from qiime.util import make_option

from cogent.parse.tree import DndParser

from qiime.pycogent_backports.parse_biom import parse_biom_table
from qiime.util import parse_command_line_parameters
from qiime.format import format_otu_table
import qiime.pycogent_backports.rich_otu_table as rich_otu_table

from qiime.simsam import sim_otu_table

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
    otu_table = parse_biom_table(otu_table_fh)
    tree_fh = open(opts.tree_file,'U')
    tree = DndParser(tree_fh)

    res_sam_names, res_otus, res_otu_mtx, res_otu_metadata = \
     sim_otu_table(otu_table.SampleIds, otu_table.ObservationIds, otu_table.iterSamples(), 
                   otu_table.ObservationMetadata, tree, opts.num, opts.dissim)


    rich_table = rich_otu_table.table_factory(res_otu_mtx,res_sam_names,res_otus,
    observation_metadata=res_otu_metadata)
    out_fh.write(rich_table.getBiomFormatJsonString())


if __name__ == "__main__":
    main()