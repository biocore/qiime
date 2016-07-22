#!/usr/bin/env python

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Will Van Treuren, Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"

from numpy import asarray
from qiime.parse import parse_mapping_file_to_dict
from collections import defaultdict

"""
Library script for creating a bipartite network that connects OTUs to Samples
and creates files suitable for graphing in Cytoscape. The bipartite graph has
been utilized for a variety of high profile papers. Usually these graphs have
shown the connections between samples and OTUs and have been arranged with a
spring embedded or force directed layout where the edge weights correspond to
the abundance of the OTU in a given sample. The clustering patterns observed
have generally corresponded to PCoA clustering patterns.

The original code that performed this function was developed by Julia Goodrich.
It was titled make_otu_network.py and is available in qiime 1.6 and earlier. It
was replaced because changing dependencies and use cases made the old code hard
to use and significantly harder to maintain.
"""


def make_sample_node_table(bt, mf_dict):
    '''Make sample node table where samples and their metadata are recorded.
    Sample node table will take the following form:

    NodeID  NodeType    Abundance  md1  md2 ...
    sample1 sample 1000    USA Europe

    where md1 and md2 ... are the metadata categories in the mapping file dict.
    inputs are biom table and mapping file dictionary created with
    parse_mapping_file_to_dict.
    '''
    # make sure to use only the sample ids found in the biom file as the
    # mapping file may have a superset of the ids found in the biom file
    sids = bt.ids()
    header = '#NodeID\tNodeType\tAbundance\t' + \
        '\t'.join(mf_dict[sids[0]].keys())
    lines = [header] + ['%s\tsample\t%s\t' % (sid, bt.data(sid, axis='sample').sum()) +
                        '\t'.join(mf_dict[sid].values()) for sid in sids]
    return lines


def make_otu_node_table(bt, md_key, md_fields):
    '''Make otu node table where nodes and their data are recorded.
    Metadata can be coded in the biom file in 3 ways, as a string, as a list,
    or as a dictionary:
    -- as a string -- {'taxonomy':'k__bacteria;p__xyz'}
    -- as a list -- {'taxonomy':['k__bacteria', 'p__xyz']}
    -- as a dictionary -- {'taxonomy':{'kingdom':'bacteria', 'phylum':'xyz'}}
    This function will attempt to divine which type of metadata exists in the
    file. If the metadata is a string, it will be assumed that it is delimited
    by semicolons. If it is not this function will fail not split it.
    Inputs:
     bt - biom table
     md_key - str, key to access the metadata str/list/dict
     md_fields - list of strs, the fields to pull from if the
     metadata is presented as a dictionary. The md_fields will be the header
     keys in the output otu_node_table. If metadata is a string or a dict then
     the number of md_fields can differ from the number of available metadata
     entries. This can cause the
     resulting file to have jagged edges (the header is longer than the metadata
     or vice versa). This does not appear to be a problem with Cytoscape. It may
     cause an error in QIIME though so we apply the same solution as used in
     summarize_taxa.py The function that takes care of it in that case is called
     sum_counts_by_consensus and it appends 'Other' to the metadata until the
     required length is reached.
    Outputs:
     list of lines of following form:
     #NodeID    NodeType    Abundance   md_field[0] md_field[1] ...
     otu1   otu 23  x   y
    '''
    header = '#NodeID\tNodeType\tAbundance\t' + '\t'.join(md_fields)
    lines = [header]
    # assume that all metadata has same format so testing any entry sufficient
    md_type = type(bt.metadata(axis='observation')[0][md_key])
    if md_type is str or md_type is unicode:
        # there are a huge number of possible ways in which the string could be
        # formatted. if its not splittable on a semicolon (preferred for qiime)
        # no splitting will occur.
        for i, otu in enumerate(bt.ids(axis='observation')):
            line = '%s\totu\t%s\t' % (otu, bt.data(otu, 'observation').sum())
            line += bt.metadata(axis='observation')[i][md_key].replace(';',
                                                                       '\t')
            lines.append(line)
    if md_type is list:
        for i, otu in enumerate(bt.ids(axis='observation')):
            line = '%s\totu\t%s\t' % (otu, bt.data(otu, 'observation').sum())
            line += '\t'.join(bt.metadata(axis='observation')[i][md_key])
            lines.append(line)
    if md_type is defaultdict:
        # if md_type is defaultdict keys in md_fields that fail will produce
        # empty lists or strs. these will cause TypeErrors in join.
        try:
            for i, otu in enumerate(bt.ids(axis='observation')):
                line = '%s\totu\t%s\t' % (otu, bt.data(otu, 'observation').sum())
                line += '\t'.join(
                    [bt.metadata(axis='observation')[i][md_key][k]
                     for k in md_fields])
                lines.append(line)
        except TypeError:
            raise ValueError('The md_fields provided were not all found in ' +
                             'the input biom table metada.')
    if md_type is dict:
        # md_fields not found will cause keyerrors
        try:
            for i, otu in enumerate(bt.ids(axis='observation')):
                line = ('%s\totu\t%s\t'
                        % (otu, bt.data(otu, 'observation').sum()))
                line += '\t'.join(
                    [bt.metadata(axis='observation')[i][md_key][k]
                     for k in md_fields])
                lines.append(line)
        except KeyError:
            raise ValueError('The md_fields provided were not all found in ' +
                             'the input biom table metada.')

    # given that the md was list or string we need to check if jagged edges
    # would be returned. instead of erroring if the output would have jagged
    # edges we apply the solution found in summarize_taxa.py. we append 'Other'
    # to the string  until it is of desired length. if the length of the
    # metadata is longer than the header we will discard \t separated entries
    # until we reach the dlen.
    if md_type is list or md_type is str or md_type is unicode:
        dlen = len(header.split('\t'))
        checked_lines = [header]
        for line in lines[1:]:  # skip header, we know its the right length
            sline = line.split('\t')
            diff = dlen - len(sline)
            if diff == 0:
                checked_lines.append(line)
            elif diff > 0:  # header longer than line, must add 'Other'*diff
                checked_lines.append(line + '\tOther' * diff)
            elif diff < 0:  # header shorter than line, remove extra entries
                checked_lines.append('\t'.join(sline[:dlen]))
        lines = checked_lines
    return lines


def make_node_attr_table(otu_node_lines, sample_node_lines,
                         sample_color, otu_color, sample_size, otu_size, sample_shape,
                         otu_shape):
    '''Make a preference table to load as node attributes for cytoscape.
    This file makes it easy to color, shape, and size the nodes according
    to the desire of the user. The color, size, and shape inputs are lists
    of strings that specify fields in the headers of otu_node_lines
    and sample_node_lines. The output will be as follows:

    #NodeID NodeType    Abundance   Color   Shape   Size
    otu1    otu 45  Bacteria_bacteriodales    spc56   xyz
    sample1 sample  56  post_treatment  tp_5    abc

    In the above example the user has passed sample_color as ['Treatment']
    and sample1 happens to be post treatment. For otu_color they passed
    ['kingdom', 'phylum'] and otu1 had kingdom Bacteria and phylum
    bacteriodales. This allows arbitrary numbers of color, size, shape
    combos to be created so that everything is fully customizable. If more
    than one field is passed the values for those fields will be joined
    with a '_'.
    Inputs:
     otu_node_lines - list of strs, output of make_otu_node_table
     sample_node_lines - list of strs, output of make_sample_node_table
     _colors, _size, _shape - each of these 6 fields must be a list of
     strings which identify which header fields are desired.
    '''
    # no comments
    sample_nodes = parse_mapping_file_to_dict(sample_node_lines)[0]
    otu_nodes = parse_mapping_file_to_dict(otu_node_lines)[0]
    header = '#NodeID\tNodeType\tAbundance\tColor\tSize\tShape'
    lines = [header]
    # make list of nodes that includes samples and otus
    nodes = sample_nodes.keys() + otu_nodes.keys()
    # make 5 lists which will be the columns of the output file
    nids, nodetypes, abundances, colors, sizes, shapes = [], [], [], [], [], []
    for node in nodes:
        if node in otu_nodes:
            nodetype_val = 'otu'
            abundance_val = otu_nodes[node]['Abundance']
            color_val = '_'.join([otu_nodes[node][i] for i in otu_color])
            size_val = '_'.join([otu_nodes[node][i] for i in otu_size])
            shape_val = '_'.join([otu_nodes[node][i] for i in otu_shape])
        elif node in sample_nodes:
            nodetype_val = 'sample'
            abundance_val = sample_nodes[node]['Abundance']
            color_val = '_'.join([sample_nodes[node][i] for i in sample_color])
            size_val = '_'.join([sample_nodes[node][i] for i in sample_size])
            shape_val = '_'.join([sample_nodes[node][i] for i in sample_shape])
        nids.append(node)
        nodetypes.append(nodetype_val)
        abundances.append(abundance_val)
        colors.append(color_val)
        sizes.append(size_val)
        shapes.append(shape_val)
    nls = ['\t'.join(vals) for vals in zip(nids, nodetypes, abundances, colors,
                                           sizes, shapes)]
    return lines + nls


def make_edge_table(bt):
    '''Make edge table where each sample is connected to the otus found in it.
    The edge table will take the following form:

    #Sample  OTU Abundance
    sample1 otu1    35

    The abundance is occurrence of the OTU and will be used to weight the edges.
    Input is a biom table.
    '''
    data = asarray([d for d in bt.iter_data(axis='observation', dense=True)])
    oids = asarray(bt.ids(axis='observation'))
    header = '#Sample\tOTU\tAbundance'
    lines = [header]
    for sample in bt.ids():
        sample_ind = bt.index(sample, 'sample')
        otu_ids = oids[data[:, sample_ind].nonzero()[0]]
        otu_abs = data[:, sample_ind][data[:, sample_ind].nonzero()[0]]
        connections = ['%s\t%s\t%s' % (sample, otu, ab) for otu, ab in
                       zip(otu_ids, otu_abs)]
        lines.extend(connections)
    return lines


def _write_table(lstrs, fp):
    '''Write a table. Input is list of strings and a filepath. Untested.'''
    o = open(fp, 'w')
    o.writelines('\n'.join(lstrs))
    o.close()
