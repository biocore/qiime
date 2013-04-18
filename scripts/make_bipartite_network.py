#!/usr/bin/env python
# File created on 17 Apr 2013
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"


from qiime.util import parse_command_line_parameters, make_option, create_dir
import os
from qiime.parse import parse_mapping_file_to_dict
from biom.parse import parse_biom_table
from qiime.make_bipartite_network import (make_sample_node_table, _write_table,
    make_otu_node_table, make_node_attr_table, make_edge_table)

script_info = {}
script_info['brief_description'] = """This script makes a bipartite network connecting samples to OTUs. It is most suitable for visualization with cytoscape."""
script_info['script_description'] = \
"""This script was created to ease the process of making bipartite networks that have appeared in high profile publications including 10.1073/pnas.1217767110, 10.1126/science.1198719, and 10.1126/science.1198719. The script will take an biom table and a mapping file and produce an edge table which connects each sample in the biom table to the OTUs found in that sample. It is bipartite because there are two distinct node classes -- OTU and Sample. The edges are weighted by the abundance of the OTU in the sample to which it is connected. The output files of this script are intended to be loaded into Cytoscape. The EdgeTable should be uploaded first, and then the NodeAttrTable file can be uploaded as node attributes to control coloring, sizing, and shaping as the user desires. The overall idea behind this script is to make bipartite network creation easier. To that end, the color, size, and shape options are used to provide fields in the NetworkViz tab of Cytoscape so that nodes can be appropriately presented. Those options are passed via comma separated strings (as in the example below). The most common visualization strategy is to color sample nodes by a metadata category like timepoint or pH, color OTU nodes by one of their taxonomic levels, and to scale OTU node size by abundance. This script makes this process easy (as well as a myriad of other visualiation strategies). Once the tables are created by this script they must be opened in cytoscape. This process is described in detail in the QIIME bipartite network tutorial available at: http://qiime.org/tutorials/making_cytoscape_networks.html All color, size, and shape options in this script default to 'NodeType'. OTU nodes have NodeType: otu, sample nodes have NodeType: sample. Thus, if you ran this script with defaults, you would only be able to change the shape, size, and color of the nodes depending on whether or not they were OTUs or samples. You would not be able to distinguish between two OTUs based on color, shape, or size. The script is flexible in that it allows you to pass any number of fields for the --{s,o}{shape,size,color}. This will allow you to distinguish between OTU and sample nodes in a huge number of different ways. The usage examples below show some of the common use cases and what options you would pass to emulate them. There are a couple of important considerations for using this script:

Note that the --md_fields option is almost required as it determines what OTU metadata will go into the OTU headers. If it is not passed the script will attempt to coerce the metadata into 7 fields (kingdom, phylum, class, ...). The script supports metadata strings, lists, dicts, and defaultdicts. With strings the levels must be separated by semicolons. The script will error if a jagged output would be produced, i.e. the number of headers would be unequal to the number of values for the metadata. In cases where the metadata is a string or list, all that is required is that the md_fields input be of equal length and in the order of the metadata that is desired.
The available fields for both sample and otu nodes are:
[NodeType, Abundance]

For otu nodes the additional fields available are:
[kingdom, phylum, class, order, family, genus, species]

For sample nodes the additional fields available are any fields found in the mapping file headers. 

If multiple fields are passed for a given option, they will be concatenated in the output with a '_' character. 
"""
script_info['script_usage'] = [\
("","""Create an EdgeTable and NodeAttrTable that allow you to color sample nodes with one of their metadata categories (Treatment for our example), OTU nodes by their taxonomic level (class for our example), control OTU node size by their abundance, and control node shape by whether its an OTU or sample.""",
"""%prog -i otu_table.biom -m mapping_file.txt -t taxonomy --md_fields 'k,p,c,o,f' -o bipartite_network/ --scolors 'Treatment' --ocolors 'c' --osize 'Abundance'"""),
("","""Create an EdgeTable and NodeAttrTable that allow you to color sample nodes by a combination of their time point and diet, color OTU nodes by their abundance and family, and node shape by whether the node is an OTU or sample. Note that the names in the --md_fields are irrelevant as long as the field passed for --ocolors is available. The length is important however, since there are 5 levels in our OTU table. If fewer fewer than 5 fields were passed for --md_fields we would get an error.""",
"""%prog -i otu_table.biom -m mapping_file.txt -t Taxonomy --md_fields 'a1,a2,a3,a4,f'-o bipartite_network/ --scolors 'TimePt,Diet' --ocolors 'f,Abundance'""")]
script_info['output_description']= """The output of this script is four files:
    1. EdgeTable - table with connections between samples and otus. 
    2. OTUNodeTable - table with OTUs and their associated metadata. 
    3. SampleNodeTable - table with samples and their associated metadata.
    4. NodeAttrTable - table with the node attributes specified by the user with
       the given options."""
script_info['required_options'] = [\
 make_option('-i','--biom_fp',type="existing_filepath",help='the input file ' +\
    'path for biom table.'),
 make_option('-m','--map_fp',type="existing_filepath",help='the input file ' +\
    'path for mapping file.'),
 make_option('-o','--out_fp',type="new_filepath",help='directory that ' +\
    'will be created for storing the results.'),
 make_option('-t','--metadata_string',type="string",help='Key to retrieve ' +\
    'taxonomy from the biom file.')]
script_info['optional_options'] = [\
 make_option('--scolors',type="string",help='commas separated string ' +\
    'specifying fields of interest for sample node coloring [default: %default].',
    default='NodeType'),
 make_option('--ocolors',type="string",help='commas separated string ' +\
    'specifying fields of interest for otu node coloring [default: %default].',
    default='NodeType'),
 make_option('--sshapes',type="string",help='commas separated string ' +\
    'specifying fields of interest for sample node shape [default: %default].',
    default='NodeType'),
 make_option('--oshapes',type="string",help='commas separated string ' +\
    'specifying fields of interest for otu node shape [default: %default].',
    default='NodeType'),
 make_option('--ssizes',type="string",help='commas separated string ' +\
    'specifying fields of interest for sample node size [default: %default].',
    default='NodeType'),
 make_option('--osizes',type="string",help='commas separated string ' +\
    'specifying fields of interest for otu node size [default: %default].',
    default='NodeType'),
 make_option('--md_fields',type="string",help='metadata fields that will be '+\
    'the headers of the OTUNodeTable. In addition, if the biom table has ' +\
    'metadata dictionaries, md_fields will be the keys extracted from the ' +\
    'biom table metadata [default: %default].',
    default='kingdom,phylum,class,order,family,genus,species')]

script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    otu_table_fp = opts.biom_fp
    map_fp = opts.map_fp
    out_fp = opts.out_fp
    taxonomy_key = opts.metadata_string
    scolors = opts.scolors.split(',')
    ocolors = opts.ocolors.split(',')
    sshapes = opts.sshapes.split(',')
    oshapes = opts.oshapes.split(',')
    ssizes = opts.ssizes.split(',')
    osizes = opts.osizes.split(',')
    md_fields = opts.md_fields.split(',')

    # check that the otu fields asked for are available
    shared_options = ['NodeType','Abundance']
    if not all([i in md_fields+shared_options for i in ocolors+oshapes+osizes]):
        raise ValueError('The fields specified for otu colors, sizes, or '+\
            'shapes are not in either the shared options (NodeType,Abundance)'+\
            ' or the supplied md_fields. These fields must be a subset of '+\
            'the union of these sets. Have you passed ocolors, osizes or '+\
            'oshapes that are not in the md_fields?')
    # check that the sample fields asked for are available. mapping file 
    # elements should all have same metadata keys
    sopts = parse_mapping_file_to_dict(map_fp)[0].items()[0][1].keys()
    if not all([i in sopts+shared_options for i in scolors+sshapes+ssizes]):
        raise ValueError('The fields specified for sample colors, sizes, or '+\
            'shapes are not in either the shared options (NodeType,Abundance)'+\
            ' or the supplied mapping file. These fields must be a subset of '+\
            'the union of these sets. Have you passed scolors, ssizes or '+\
            'sshapes that are not in the mapping file headers?')

    # actual compuation begins
    create_dir(out_fp, fail_on_exist=True)

    bt = parse_biom_table(open(otu_table_fp))
    pmf = parse_mapping_file_to_dict(map_fp)[0] # [1] is comments, don't need
    sample_node_table = make_sample_node_table(bt, pmf)
    otu_node_table = make_otu_node_table(bt, opts.metadata_string, md_fields)
    node_attr_table = make_node_attr_table(otu_node_table, sample_node_table,
        scolors, ocolors, ssizes, osizes, sshapes, oshapes)
    edge_table = make_edge_table(bt)

    _write_table(sample_node_table, os.path.join(out_fp,'SampleNodeTable.txt'))
    _write_table(otu_node_table, os.path.join(out_fp,'OTUNodeTable.txt'))
    _write_table(node_attr_table, os.path.join(out_fp,'NodeAttrTable.txt'))
    _write_table(edge_table, os.path.join(out_fp,'EdgeTable.txt'))

if __name__ == "__main__":
    main()