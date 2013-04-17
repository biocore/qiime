#!/usr/bin/env python

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2013, The QIIME Project" 
__credits__ = ["Will Van Treuren, Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

from numpy import array
from qiime.parse import parse_mapping_file_to_dict

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
    sids = bt.SampleIds
    header = '#NodeID\tNodeType\tAbundance\t'+'\t'.join(mf_dict[sids[0]].keys())
    lines = [header] + ['%s\tsample\t%s\t' % (sid, bt.sampleData(sid).sum()) +
            '\t'.join(mf_dict[sid].values()) for sid in sids]
    return lines

def make_otu_node_table(bt, md_key, num_lvls=7):
    '''Make otu node table where nodes and their data are recorded.
    Otu node table will take the following form:
    
    NodeID  NodeType    Abundance  kingdom  phylum ...
    otu1    otu    35  bacteria    bacteroides

    Inputs:
     bt - biom table
     md_key - str, key to retrieve taxonomy info from the biom table
     lvls - int, default 7, number of levels in the taxonomy to record. useful
     if the table only has taxonomy to x<7 levels.
    ''' 
    lvls = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 
        'species'][:num_lvls]
    header = '#NodeID\tNodeType\tAbundance\t'+'\t'.join(lvls)
    lines = [header]

    # metadata can be of a variety of types: 
    # it can be a single string - {'taxonomy':'k__bacteria;p__xyz'} 
    # a list - {'taxonomy':['k__bacteria', 'p__xyz']}
    # or a dict of values - {'taxonomy':{'kingdom':'bacteria', 'phylum':'xyz'}}
    # we will only support a dict or list. with dict the keys must be specific
    if type(bt.ObservationMetadata[0][md_key]) == dict:
        # require that the keys for the dict be kingdom, phylum, class ...
        # otherwise have to support insane number of different specifications
        if all([i in lvls for i in bt.ObservationMetadata[0][md_key].keys()]):
            for i,otu in enumerate(bt.ObservationIds):
                l1 = '%s\totu\t%s\t' % (otu, bt.observationData(otu).sum()) 
                l2 = '\t'.join([bt.ObservationMetadata[i][md_key][lvl] for lvl 
                    in lvls])
                lines.append(l1+l2)
        else:
            raise ValueError('The biom table has observation metadata that '+\
                'keyed on something that is not kingdom, phylum, class... '+\
                'This library only supports metadata keyed this way. For '+\
                'example: taxonomy:{kingdom:bacteria,phylum:something,...}.')
    elif type(bt.ObservationMetadata[0][md_key]) == list:
        for i,otu in enumerate(bt.ObservationIds):
            line = '%s\totu\t%s\t' % (otu, bt.observationData(otu).sum()) + \
                '\t'.join(bt.ObservationMetadata[i][md_key])
            lines.append(line)
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
    with a _ character (as in the otu example above). 
    Inputs:
     otu_node_lines - list of strs, output of make_otu_node_table
     sample_node_lines - list of strs, output of make_sample_node_table
     _colors, _size, _shape - each of these 6 fields must be a list of 
     strings which identify which header fields are desired. 
    '''
    sample_nodes = parse_mapping_file_to_dict(sample_node_lines)[0] #no comments
    otu_nodes = parse_mapping_file_to_dict(otu_node_lines)[0] 
    #print sample_nodes, otu_node
    header = '#NodeID\tNodeType\tAbundance\tColor\tSize\tShape'
    lines = [header]
    # make list of nodes that includes samples and otus
    nodes = sample_nodes.keys()+otu_nodes.keys()
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
    return lines+nls

def make_edge_table(bt):
    '''Make edge table where each sample is connected to the otus found in it.
    The edge table will take the following form:
    
    #Sample  OTU Abundance   
    sample1 otu1    35

    The abundance is occurrence of the OTU and will be used to weight the edges.
    Input is a biom table.
    '''
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    oids = array(bt.ObservationIds)
    header = '#Sample\tOTU\tAbundance'
    lines = [header]
    for sample in bt.SampleIds:
        sample_ind = bt.getSampleIndex(sample)
        otu_ids = oids[data[:,sample_ind].nonzero()[0]]
        otu_abs = data[:,sample_ind][data[:,sample_ind].nonzero()[0]]
        connections = ['%s\t%s\t%s' % (sample, otu, ab) for otu,ab in 
            zip(otu_ids, otu_abs)]
        lines.extend(connections)
    return lines

def _write_table(lstrs, fp):
    '''Write a table. Input is list of strings and a filepath.'''
    o = open(fp, 'w')
    o.writelines('\n'.join(lstrs))
    o.close()