#!/usr/bin/env python
# File created on 17 Apr 2013
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren, Julia Goodrich, Luke Ursell"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"


import os

from biom import load_table

from qiime.util import parse_command_line_parameters, make_option, create_dir
from qiime.parse import parse_mapping_file_to_dict
from qiime.make_bipartite_network import (make_sample_node_table, _write_table,
                                          make_otu_node_table, make_node_attr_table, make_edge_table)

script_info = {}
script_info[
    'brief_description'] = ("This script makes a bipartite network connecting "
                            "samples to observations. It is most suitable for "
                            "visualization with cytoscape.")
script_info['script_description'] = \
    """This script was created to ease the process of making bipartite networks
that have appeared in high profile publications including
10.1073/pnas.1217767110 and 10.1126/science.1198719. The script will take a
biom table and a mapping file and produce an edge table which connects each
sample in the biom table to the observations found in that sample. It is
bipartite because there are two distinct node classes -- OTU and Sample. The
'OTU' node class does not have to be an operational taxonomic unit, it can be a
KEGG category or metabolite etc.  -- anything that is an observation. The edges
are weighted by the abundance of the observation in the sample to which it is
connected. The output files of this script are intended to be loaded into
Cytoscape. The EdgeTable should be uploaded first, and then the NodeAttrTable
file can be uploaded as node attributes to control coloring, sizing, and
shaping as the user desires. The overall idea behind this script is to make
bipartite network creation easier.  To that end, the color, size, and shape
options are used to provide fields in the NetworkViz tab of Cytoscape so that
nodes can be appropriately presented.  Those options are passed via comma
separated strings (as in the example below).  The most common visualization
strategy is to color sample nodes by a metadata category like timepoint or pH,
color OTU nodes by one of their taxonomic levels, and to scale OTU node size by
abundance. This script makes this process easy (as well as a myriad of other
visualiation strategies). Once the tables are created by this script they must
be opened in cytoscape. This process is described in detail in the QIIME
bipartite network tutorial available at:
http://qiime.org/tutorials/making_cytoscape_networks.html All color, size, and
shape options in this script default to 'NodeType'. OTU nodes have NodeType:
otu, sample nodes have NodeType: sample. Thus, if you ran this script with
defaults, you would only be able to change the shape, size, and color of the
nodes depending on whether or not they were observations or samples. You would
not be able to distinguish between two observations based on color, shape, or
size. The script is flexible in that it allows you to pass any number of fields
for the --{s,o}{shape,size,color}. This will allow you to distinguish between
OTU and sample nodes in a huge number of different ways. The usage examples
below show some of the common use cases and what options you would pass to
emulate them. There are a couple of important considerations for using this
script:

Note that the --md_fields option has a different meaning depending on the type
of metadata in the biom table. Regardless of type, the md_fields will be the
headers in the OTUNodeTable.txt. If the metadata is a dict or default dict, the
md_fields will be used as keys to extract data from the biom file metadata. If
the metadata is a list or a string, then the md_fields will be have no
intrinsic relation to the columns they head. For example if
md_fields=['k','p','c'] and the metadata contained in a given OTU was
'k__Bacteria;p__Actinobacter;c__Actino' the resulting OTUNodeTable would have
k__Bacteria in the 'k' column, p__Actinobacter in the 'p' column, and c__Actino
in the 'c' column. If one passed md_fields=['1.0','XYZ','Five'] then the
OTUNodeTable would have columns headed by ['1.0','XYZ','Five'], but the
metadata values in those columns would be the same (e.g. '1.0' column entry
would be k__Bacteria etc.) If the number of elements in the metadata for a
given OTU is not equal to the number of headers provided the script will adjust
the OTU metadata. In the case where the metadata is too short, it will add
'Other' into the OTU metadata until the required length is reached.  In the
case where the metadata is too long it will simply remove extra entries. This
means you can end up with many observations which have the value of 'Other' if
you have short taxonomic strings/lists for your observations.

The available fields for both sample and otu nodes are:
[NodeType, Abundance]

For observation nodes the additional fields available are:
any fields you passed for the md_fields

For sample nodes the additional fields available are
any fields found in the mapping file headers

If multiple fields are passed for a given option, they will be concatenated in the output with a '_' character.
"""
script_info['script_usage'] = [
    ("", ("Create an EdgeTable and NodeAttrTable that allow you to color "
          "sample nodes with one of their metadata categories (Treatment for "
          "our example), observation nodes (in this case OTUs) by their "
          "taxonomic level (class for our example), control observation node "
          "size by their abundance, and control node shape by whether its an "
          "observation or sample."),
    ("%prog -i otu_table.biom -m mapping_file.txt -k taxonomy --md_fields "
     "'k,p,c,o,f' -o bipartite_network/ --scolors 'Treatment' --ocolors 'c' "
     "--osize 'Abundance'")),
    ("", ("Create an EdgeTable and NodeAttrTable that allow you to color "
          "sample nodes by a combination of their time point and diet, color "
          "observation nodes by their abundance and family, and node shape by "
          "whether the node is an observation or sample. Note that the names "
          "in the --md_fields are irrelevant as long as the field passed for "
          "--ocolors is available. The length is important however, since "
          "there are 5 levels in our OTU table. If fewer fewer than 5 fields "
          "were passed for --md_fields we would get an error."),
    ("%prog -i otu_table.biom -m mapping_file.txt -k taxonomy --md_fields "
     "'a1,a2,a3,a4,f' -o bipartite_network_combo_colors/ --scolors "
     "'TimePt,Diet' --ocolors 'f,Abundance'"))]
script_info[
    'script_usage_output_to_remove'] = [
    '$PWD/bipartite_network/',
    '$PWD/bipartite_network_combo_colors/']
script_info['output_description'] = """The output of this script is four files:
    1. EdgeTable - table with connections between samples and observations.
    2. OTUNodeTable - table with observations and their associated metadata.
    3. SampleNodeTable - table with samples and their associated metadata.
    4. NodeAttrTable - table with the node attributes specified by the user with
       the given options."""
script_info['required_options'] = [
    make_option('-i', '--biom_fp', type="existing_filepath",
                help='the input file path for biom table.'),
    make_option('-m', '--map_fp', type="existing_filepath",
                help='the input file path for mapping file.'),
    make_option('-o', '--output_dir', type="new_dirpath",
                help='directory to be created for storing the results.'),
    make_option('-k', '--observation_md_header_key', type="string",
                help=("Key to retrieve metadata (usually taxonomy) from the "
                      "biom file.")),
    make_option('--md_fields', type="string", help=('metadata fields that '
                'will be the headers of the OTUNodeTable. If the biom table '
                'has metadata dictionaries, md_fields will be the keys '
                'extracted from the biom table metadata. Passed like '
                '"kingdom,phylum,class".'))]
script_info['optional_options'] = [
    make_option('--scolors', type="string", help=('commas separated string '
                'specifying fields of interest for sample node coloring '
                '[default: %default].'),
                default='NodeType'),
    make_option('--ocolors', type="string", help=('commas separated string '
                'specifying fields of interest for observation node coloring '
                '[default: %default].'),
                default='NodeType'),
    make_option('--sshapes', type="string", help=('commas separated string '
                'specifying fields of interest for sample node shape '
                '[default: %default].'),
                default='NodeType'),
    make_option('--oshapes', type="string", help=('commas separated string '
                'specifying fields of interest for observation node shape '
                '[default: %default].'),
                default='NodeType'),
    make_option('--ssizes', type="string", help=('commas separated string '
                'specifying fields of interest for sample node size '
                '[default: %default].'),
                default='NodeType'),
    make_option('--osizes', type="string", help=('commas separated string '
                'specifying fields of interest for observation node size '
                '[default: %default].'),
                default='NodeType')]


script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    otu_table_fp = opts.biom_fp
    map_fp = opts.map_fp
    output_dir = opts.output_dir
    scolors = opts.scolors.split(',')
    ocolors = opts.ocolors.split(',')
    sshapes = opts.sshapes.split(',')
    oshapes = opts.oshapes.split(',')
    ssizes = opts.ssizes.split(',')
    osizes = opts.osizes.split(',')
    md_fields = opts.md_fields.split(',')

    # check that the otu fields asked for are available
    shared_options = ['NodeType', 'Abundance']
    if not all([i in md_fields + shared_options for i in ocolors + oshapes + osizes]):
        option_parser.error('The fields specified for observation colors, '
                            'sizes, or shapes are not in either the shared '
                            'options (NodeType,Abundance) or the supplied '
                            'md_fields. These fields must be a subset of the '
                            'union of these sets. Have you passed ocolors, '
                            'osizes or oshapes that are not in the md_fields?')
    # check that the sample fields asked for are available. mapping file
    # elements should all have same metadata keys
    sopts = parse_mapping_file_to_dict(map_fp)[0].items()[0][1].keys()
    if not all([i in sopts + shared_options for i in scolors + sshapes + ssizes]):
        option_parser.error('The fields specified for sample colors, sizes, '
                            'or shapes are not in either the shared options '
                            '(NodeType,Abundance) or the supplied mapping '
                            'file. These fields must be a subset of the union '
                            'of these sets. Have you passed scolors, ssizes '
                            'or sshapes that are not in the mapping file '
                            'headers?')

    # actual compuation begins
    try:
        create_dir(output_dir, fail_on_exist=True)
    except OSError:
        option_parser.error('Directory already exists. Will not overwrite.')

    bt = load_table(otu_table_fp)
    pmf = parse_mapping_file_to_dict(map_fp)[0]  # [1] is comments, don't need
    sample_node_table = make_sample_node_table(bt, pmf)
    otu_node_table = make_otu_node_table(bt, opts.observation_md_header_key,
                                         md_fields)
    node_attr_table = make_node_attr_table(otu_node_table, sample_node_table,
                                           scolors, ocolors, ssizes, osizes, sshapes, oshapes)
    edge_table = make_edge_table(bt)

    _write_table(
        sample_node_table,
        os.path.join(
            output_dir,
            'SampleNodeTable.txt'))
    _write_table(otu_node_table, os.path.join(output_dir, 'OTUNodeTable.txt'))
    _write_table(
        node_attr_table,
        os.path.join(
            output_dir,
            'NodeAttrTable.txt'))
    _write_table(edge_table, os.path.join(output_dir, 'EdgeTable.txt'))

if __name__ == "__main__":
    main()
