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


from qiime.util import parse_command_line_parameters, make_option
import os
from qiime.parse import parse_mapping_file_to_dict
from biom.parse import parse_biom_table
from qiime.make_bipartite_network import (make_sample_node_table, _write_table,
    make_otu_node_table, make_node_attr_table, make_edge_table)

script_info = {}
script_info['brief_description'] = """This script makes a bipartite network 
    connecting samples to OTUs. It is most suitable for visualization with 
    cytoscape."""
script_info['script_description'] = \
"""This script was created to ease the process of making bipartite networks that have appeared in high profile publications including 10.1073/pnas.1217767110, 10.1126/science.1198719, and 10.1126/science.1198719. The script will take an biom table and a mapping file and produce an edge table which connects each sample in the biom table to the OTUs found in that sample. It is bipartite because there are two distinct node classes -- OTU and Sample. The edges are weighted by the abundance of the OTU in the sample to which it is connected. The output files of this script are intended to be loaded into Cytoscape. The EdgeTable should be uploaded first, and then the NodeAttrTable file can be uploaded as node attributes to control coloring, sizing, and shaping as the user desires. The overall idea behind this script is to make bipartite network creation easier. To that end, the color, size, and shape options are used to provide fields in the NetworkViz tab of Cytoscape so that nodes can be appropriately presented. Those options are passed via comma separated strings (as in the example below). The most common visualization strategy is to color sample nodes by a metadata category like timepoint or pH, color OTU nodes by one of their taxonomic levels, and to scale OTU node size by abundance. This script makes this process easy (as well as a myriad of other visualiation strategies). Once the tables are created by this script they must be opened in cytoscape. This process is described in detail in the QIIME bipartite network tutorial available at: http://qiime.org/tutorials/making_cytoscape_networks.html All color, size, and shape options in this script default to 'NodeType'. OTU nodes have NodeType: otu, sample nodes have NodeType: sample. Thus, if you ran this script with defaults, you would only be able to change the shape, size, and color of the nodes depending on whether or not they were OTUs or samples. You would not be able to distinguish between two OTUs based on color, shape, or size. The script is flexible in that it allows you to pass any number of fields for the --{s,o}{shape,size,color}. This will allow you to distinguish between OTU and sample nodes in a huge number of different ways. The usage examples below show some of the common use cases and what options you would pass to emulate them. There are a couple of important considerations for using this script:

The available fields for both sample and otu nodes are:
[NodeID, NodeType, Abundance]

For otu nodes the additional fields available are:
[kingdom, phylum, class, order, family, genus, species]

For sample nodes the additional fields available are any field found in the mapping file. 

If multiple fields are passed for a given option, they will be concatenated in the header with a '_' character. Their values for the corresponding sample or OTU will also be concatenated with a '_'.

If the biom table provided has dictionaries of metadata, the dictionaries must be keyed by kingdom, phylum, class ... etc. or the script will error. If the metadata is provided as a list, it should be in the order of kingdom, phylum, ... If it is not, the fields in the otu node table will be improperly written. If the metadata is provided as a string the script will attempt to split on semicolons. If it fails to do this it will error."""
script_info['script_usage'] = [\
("","""Create an EdgeTable and NodeAttrTable that allow you to color sample nodes with one of their metadata categories (pH for our example), OTU nodes by their taxonomic level (class for our example), control OTU node size by their abundance, and control node shape by whether its an OTU or sample.""",
"""%prog -i otu_table.biom -m mapping_file.txt -t Taxonomy -o bipartite_network/ --scolors 'pH' --ocolors 'class' --osize 'Abundance'"""),
("","""Create an EdgeTable and NodeAttrTable that allow you to color sample nodes by a combination of their Timepoint and pH, color OTU nodes by their abundance, and node shape by whether the node is an OTU or sample.""",
"""%prog -i otu_table.biom -m mapping_file.txt -t Taxonomy -o bipartite_network/ --scolors 'Timepoint,pH' --ocolors 'Abundance'""")]
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
 make_option('--sshape',type="string",help='commas separated string ' +\
    'specifying fields of interest for sample node shape [default: %default].',
    default='NodeType'),
 make_option('--oshape',type="string",help='commas separated string ' +\
    'specifying fields of interest for otu node shape [default: %default].',
    default='NodeType'),
 make_option('--ssize',type="string",help='commas separated string ' +\
    'specifying fields of interest for sample node size [default: %default].',
    default='NodeType'),
 make_option('--osize',type="string",help='commas separated string ' +\
    'specifying fields of interest for otu node size [default: %default].',
    default='NodeType'),
 make_option('--taxonomic_levels',type="int",help='integer specifying the number' +\
    ' of levels of taxonomy to record in the OTUNodeTable. Helpful if there are'+\
    ' fewer than 7 levels available [default: %default].',
    default=7)]

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
    sshape = opts.sshape.split(',')
    oshape = opts.oshape.split(',')
    ssize = opts.ssize.split(',')
    osize = opts.osize.split(',')
    lvl = opts.taxonomic_levels

    os.mkdir(opts.out_fp)

    bt = parse_biom_table(open(otu_table_fp))
    pmf = parse_mapping_file_to_dict(map_fp)[0] # [1] is comments, don't need
    sample_node_table = make_sample_node_table(bt, pmf)
    otu_node_table = make_otu_node_table(bt, opts.metadata_string, lvl)
    node_attr_table = make_node_attr_table(otu_node_table, sample_node_table,
        scolors, ocolors, ssize, osize, sshape, oshape)
    edge_table = make_edge_table(bt)

    _write_table(sample_node_table, out_fp+'SampleNodeTable.txt')
    _write_table(otu_node_table, out_fp+'OTUNodeTable.txt')
    _write_table(node_attr_table, out_fp+'NodeAttrTable.txt')
    _write_table(edge_table, out_fp+'EdgeTable.txt')

if __name__ == "__main__":
    main()