#!/usr/bin/env python
# File created on 04 Dec 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from biom.parse import parse_biom_table
from qiime.format import format_biom_table
from qiime.filter import get_otu_ids_from_taxonomy_f
from qiime.util import parse_command_line_parameters, make_option

script_info = {}
script_info['brief_description'] = "Filter taxa from an OTU table"
script_info['script_description'] = "This scripts filters an OTU table based on taxonomic metadata. It can be applied for positive filtering (i.e., keeping only certain taxa), negative filtering (i.e., discarding only certain taxa), or both at the same time."
script_info['script_usage'] = []
script_info['script_usage'].append(("","Filter otu_table.biom to include only OTUs identified as __Bacteroidetes or p__Firmicutes.","%prog -i otu_table.biom -o otu_table_bac_firm_only.biom -p p__Bacteroidetes,p__Firmicutes"))
script_info['script_usage'].append(("","Filter otu_table.biom to exclude OTUs identified as p__Bacteroidetes or p__Firmicutes.","%prog -i otu_table.biom -o otu_table_non_bac_firm.biom -n p__Bacteroidetes,p__Firmicutes"))
script_info['script_usage'].append(("","Filter otu_table.biom to include OTUs identified as p__Firmicutes but not c__Clostridia.","%prog -i otu_table.biom -o otu_table_all_firm_but_not_clos.biom -p p__Firmicutes -n c__Clostridia"))

script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_otu_table_fp',
             type="existing_filepath",
             help='the input otu table filepath'),
 make_option('-o','--output_otu_table_fp',
             type="new_filepath",
             help='the output otu table filepath'),
]
script_info['optional_options'] = [
 make_option('-p','--positive_taxa',
             help='comma-separated list of taxa to retain [default: None; retain all taxa]'),
 make_option('-n','--negative_taxa',
             help='comma-separated list of taxa to discard [default: None; retain all taxa]'),
 make_option('--metadata_field',default='taxonomy',
             help='observation metadata identifier to filter based on [default: %default]'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    input_table = parse_biom_table(open(opts.input_otu_table_fp,'U'))
    output_table_f = open(opts.output_otu_table_fp,'w')
    metadata_field = opts.metadata_field
    positive_taxa = opts.positive_taxa
    negative_taxa = opts.negative_taxa
    
    if positive_taxa != None:
        positive_taxa = positive_taxa.split(',')
    else:
        positive_taxa = None

    if negative_taxa != None:
        negative_taxa = negative_taxa.split(',')
    else:
        negative_taxa = None
    
    filter_fn = get_otu_ids_from_taxonomy_f(positive_taxa,negative_taxa,metadata_field)
    output_table = input_table.filterObservations(filter_fn)
    output_table_f.write(format_biom_table(output_table))
    output_table_f.close()
    


if __name__ == "__main__":
    main()