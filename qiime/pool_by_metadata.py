#!/usr/bin/env python
from __future__ import division

import numpy
from qiime.parse import parse_mapping_file, parse_otu_table
from qiime.format import format_otu_table
from qiime.filter_by_metadata import parse_metadata_state_descriptions, get_sample_ids

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Justin Kuczynski", "Catherine Lozupone"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

def pool_map(map_infile, map_outfile,
    pooled_sample_name, sample_ids_to_pool):
    """pools map file according to specified criteria."""
    map_data, map_header, map_comments = parse_mapping_file(map_infile)
    map_infile.close()
    result = _pool_map(map_data, map_header, pooled_sample_name, sample_ids_to_pool)
    map_outfile.write('\n'.join(result))

def get_category_values(map_data, map_header, metadata_column):
    """gets the values in a specified category"""
    values = []
    column_index = map_header.index(metadata_column)
    for line in map_data:
        val = line[column_index]
        if val not in values:
            values.append(val)
    return values
    
def pool_iterative(map_infile, map_outfile, otu_infile, otu_outfile, metadata_column):
    """pools all samples in the map file that have the same value in the metadata column"""
    map_data, map_header, map_comments = parse_mapping_file(map_infile)
    map_infile.close()
    otu_table = parse_otu_table(otu_infile)
    otu_infile.close()
    new_mapping = None
    new_otu_table = None
    values = get_category_values(map_data, map_header, metadata_column)
    for i in values:
        if new_mapping:
            map_data, map_header, map_comments = parse_mapping_file(new_mapping)
        if new_otu_table:
            otu_table = parse_otu_table(new_otu_table.split('\n'))
        valid_states_str = metadata_column + ':' + i
        pooled_sample_name = metadata_column + '.' + i
        valid_states = parse_metadata_state_descriptions(valid_states_str)
        sample_ids_to_pool = get_sample_ids(map_data, map_header, valid_states)
        new_mapping = _pool_map(map_data, map_header,
            pooled_sample_name, sample_ids_to_pool)
        new_otu_table = _pool_otu_table(otu_table, pooled_sample_name, sample_ids_to_pool)
    map_outfile.write('\n'.join(new_mapping))
    otu_outfile.write(new_otu_table)

def _pool_map(map_data, map_header, pooled_sample_name, sample_ids_to_pool):
    """pools the map file according to specified criteria."""
    # valid_states = parse_metadata_state_descriptions(valid_states_str)
    # sample_ids = get_sample_ids(map_data, map_header, valid_states)

    # write out the filtered mapping file
    sample_id_idx = map_header.index('SampleID')

    # separate the samples to be pooled from the rest (new_map_data)
    new_map_data = []
    pooled_map_data = []
    for sam in map_data:
        if sam[sample_id_idx] in sample_ids_to_pool:
            pooled_map_data.append(sam)
        else:
            new_map_data.append(sam)
    
    # make the new pooled sample
    newsam = ['multipleValues'] * len(map_header)

    for i in range(len(map_header)):
        pooled_vals = [sam[i] for sam in pooled_map_data]
        if len(set(pooled_vals)) == 1:
            newsam[i] = pooled_vals[0]

    newsam[sample_id_idx] = pooled_sample_name
    
    new_map_data.append(newsam)

    header_line = '#' + '\t'.join(map_header)

    return [header_line] + map('\t'.join, new_map_data)

def pool_otu_table(otu_infile, otu_outfile, 
    pooled_sample_name, sample_ids_to_pool):
    """pools otu table file according to specified criteria."""
    ## otu table
    otu_table = parse_otu_table(otu_infile)
    otu_infile.close()
    result = _pool_otu_table(otu_table, pooled_sample_name, sample_ids_to_pool)
    otu_outfile.write(result)

def _pool_otu_table(otu_table, pooled_sample_name, sample_ids_to_pool):
    """pools otu table file according to specified criteria"""
    pool_sample_idxs = []
    nonpool_sample_idxs = []
    for i in range(len(otu_table[0])): #sample ids
        if otu_table[0][i] in sample_ids_to_pool:
            pool_sample_idxs.append(i)
        else:
            nonpool_sample_idxs.append(i)
    
    new_sample_ids = []
    for i in range(len(otu_table[0])): #sample ids
        if otu_table[0][i] not in sample_ids_to_pool: 
            # from valid_states string on mapfile
            new_sample_ids.append(otu_table[0][i])
    new_sample_ids.append(pooled_sample_name)
    
    # otu mtx
    new_sample_abund = otu_table[2][:,pool_sample_idxs].sum(1)
    newdims = (len(otu_table[2]),len(new_sample_ids))

    new_otu_mtx = numpy.zeros(newdims,dtype=otu_table[2].dtype)
    new_otu_mtx[:,:-1] = otu_table[2][:,nonpool_sample_idxs]
    new_otu_mtx[:,-1] = new_sample_abund
    
    return format_otu_table(new_sample_ids, otu_table[1], 
        new_otu_mtx, taxonomy=otu_table[3])
    
