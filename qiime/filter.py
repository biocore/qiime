#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Will Van Treuren", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from random import shuffle
from numpy import array, inf
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import parse_otu_table, parse_distmat, parse_mapping_file, parse_metadata_state_descriptions
from qiime.format import format_otu_table, format_distance_matrix, format_mapping_file

def sample_ids_from_metadata_description(mapping_f,valid_states_str):
    """ Given a description of metadata, return the corresponding sample ids
    """
    map_data, map_header, map_comments = parse_mapping_file(mapping_f)
    valid_states = parse_metadata_state_descriptions(valid_states_str)
    sample_ids = get_sample_ids(map_data, map_header, valid_states)
    return sample_ids

def get_sample_ids(map_data, map_header, states):
    """Takes col states in {col:[vals]} format.

    If val starts with !, exclude rather than include.
    
    Combines cols with and, states with or.

    For example, Study:Dog,Hand will return rows where Study is Dog or Hand;
    Study:Dog,Hand;BodySite:Palm,Stool will return rows where Study is Dog
    or Hand _and_ BodySite is Palm or Stool; Study:*,!Dog;BodySite:*,!Stool
    will return all rows except the ones where the Study is Dog or the BodySite
    is Stool.
    """
    
    name_to_col = dict([(s,map_header.index(s)) for s in states])
    good_ids = []
    for row in map_data:    #remember to exclude header
        include = True
        for s, vals in states.items():
            curr_state = row[name_to_col[s]]
            include = include and (curr_state in vals or '*' in vals) \
                and not '!'+curr_state in vals
        if include:        
            good_ids.append(row[0])
    return good_ids


def filter_fasta(input_seqs,output_seqs_f,seqs_to_keep,negate=False):
    """ Write filtered input_seqs to output_seqs_f which contains only seqs_to_keep
    """
    seqs_to_keep_lookup = {}.fromkeys([seq_id.split()[0]
                               for seq_id in seqs_to_keep])
    # Define a function based on the value of negate
    if not negate:
        def keep_seq(seq_id):
            return seq_id.split()[0] in seqs_to_keep_lookup
    else:
        def keep_seq(seq_id):
            return seq_id.split()[0] not in seqs_to_keep_lookup
    
    for seq_id, seq in input_seqs:
        if keep_seq(seq_id):
            output_seqs_f.write('>%s\n%s\n' % (seq_id, seq))
    output_seqs_f.close()

def filter_mapping_file(mapping_f,sample_ids_to_keep):
    mapping_data, header, comments = parse_mapping_file(mapping_f)
    filtered_mapping_data = []
    for entry in mapping_data:
        if entry[0] in {}.fromkeys(sample_ids_to_keep):
            filtered_mapping_data.append(entry)
    return format_mapping_file(header,filtered_mapping_data)

def filter_samples_from_distance_matrix(dm,samples_to_discard,negate=False):
    """ Remove specified samples from distance matrix 
    
        dm: (sample_ids, dm_data) tuple, as returned from 
         qiime.parse.parse_distmat; or a file handle that can be passed
         to qiime.parse.parse_distmat
    
    """
    try:
        sample_ids, dm_data = dm
    except ValueError:
        # input was provide as a file handle
        sample_ids, dm_data = parse_distmat(dm)
    
    sample_lookup = {}.fromkeys([e.split()[0] for e in samples_to_discard])
    temp_dm_data = []
    new_dm_data = []
    new_sample_ids = []
    
    if negate:
        def keep_sample(s):
            return s in sample_lookup
    else:
        def keep_sample(s):
            return s not in sample_lookup
            
    for row,sample_id in zip(dm_data,sample_ids):
        if keep_sample(sample_id):
            temp_dm_data.append(row)
            new_sample_ids.append(sample_id)
    temp_dm_data = array(temp_dm_data).transpose()
    
    for col,sample_id in zip(temp_dm_data,sample_ids):
        if keep_sample(sample_id):
            new_dm_data.append(col)
    new_dm_data = array(new_dm_data).transpose()
    
    return format_distance_matrix(new_sample_ids, new_dm_data)

def negate_tips_to_keep(tips_to_keep, tree):
    """ Return the list of tips in the tree that are not in tips_to_keep"""
    tips_to_keep = set(tips_to_keep)
    tips = set([tip.Name for tip in tree.tips()])
    return tips - tips_to_keep

def get_seqs_to_keep_lookup_from_seq_id_file(id_to_keep_f):
    """generate a lookup dict of chimeras in chimera file."""
    return set([l.split()[0].strip() for l in id_to_keep_f if not l.startswith('#') and l])
get_seq_ids_from_seq_id_file = get_seqs_to_keep_lookup_from_seq_id_file

def get_seqs_to_keep_lookup_from_fasta_file(fasta_f):
    """return the sequence ids within the fasta file"""
    return set([seq_id.split()[0] for seq_id,seq in MinimalFastaParser(fasta_f)])
get_seq_ids_from_fasta_file = get_seqs_to_keep_lookup_from_fasta_file

# start functions used by filter_samples_from_otu_table.py and filter_otus_from_otu_table.py

def get_filter_function(ids_to_keep,min_count,max_count):
    
    def f(data_vector, id_, metadata):
        return (id_ in ids_to_keep) and \
               (min_count <= data_vector.sum() <= max_count)
    return f

def filter_samples_from_otu_table(otu_table,ids_to_keep,min_count,max_count):
    filter_f = get_filter_function({}.fromkeys(ids_to_keep),
                                           min_count,
                                           max_count)
    return otu_table.filterSamples(filter_f)

def filter_otus_from_otu_table(otu_table,ids_to_keep,min_count,max_count):
    filter_f = get_filter_function({}.fromkeys(ids_to_keep),
                                           min_count,
                                           max_count)
    return otu_table.filterObservations(filter_f)

# end functions used by filter_samples_from_otu_table.py and filter_otus_from_otu_table.py



