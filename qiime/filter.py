#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from random import shuffle
from numpy import array
from qiime.parse import parse_otu_table, parse_distmat
from qiime.format import format_otu_table, format_distance_matrix

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
    
def filter_otus_from_otu_table(otu_table_lines,otus_to_discard,negate=False):
    """ Remove specified OTUs from otu_table """
    otu_table_data = parse_otu_table(otu_table_lines)
    
    otu_lookup = {}.fromkeys([e.split()[0] for e in otus_to_discard])
    new_otu_table_data = []
    new_otu_ids = []
    new_taxa = []
    
    if negate:
        def keep_otu(s):
            return s in otu_lookup
    else:
        def keep_otu(s):
            return s not in otu_lookup
    
    sample_ids, otu_ids, otu_table_data, taxa = otu_table_data
    
    for row,otu_id,taxonomy in zip(otu_table_data,otu_ids,taxa):
        if keep_otu(otu_id):
            new_otu_table_data.append(row)
            new_otu_ids.append(otu_id)
            new_taxa.append(taxonomy)
    
    new_otu_table_data = array(new_otu_table_data)
            
    result = format_otu_table(sample_ids,
                              new_otu_ids,
                              new_otu_table_data,
                              new_taxa).split('\n')
    return result
    
def filter_samples_from_otu_table(otu_table_lines,
                                  samples_to_discard,
                                  negate=False):
    """ Remove specified samples from OTU table """
    otu_table_data = parse_otu_table(otu_table_lines)
    
    sample_lookup = {}.fromkeys([e.split()[0] for e in samples_to_discard])
    new_otu_table_data = []
    new_sample_ids = []
    
    if negate:
        def keep_sample(s):
            return s in sample_lookup
    else:
        def keep_sample(s):
            return s not in sample_lookup
    
    sample_ids, otu_ids, otu_table_data, taxa = otu_table_data
    otu_table_data = otu_table_data.transpose()
    
    for row,sample_id in zip(otu_table_data,sample_ids):
        if keep_sample(sample_id):
            new_otu_table_data.append(row)
            new_sample_ids.append(sample_id)
    
    new_otu_table_data = array(new_otu_table_data).transpose()
    
    result = format_otu_table(new_sample_ids,
                              otu_ids,
                              new_otu_table_data,
                              taxa,
                              skip_empty=True).split('\n')
    return result
    
def filter_otu_table_to_n_samples(otu_table_lines,n):
    """
        randomly select n samples from the otu table
    """
    if n < 1:
        raise ValueError,\
         "number of randomly selected sample ids must be greater than 1"
    sample_ids, otu_ids, otu_table_data, taxa = parse_otu_table(otu_table_lines)
    
    samples_to_keep = list(sample_ids)
    shuffle(samples_to_keep)
    samples_to_keep = samples_to_keep[:n]
    
    otu_table_lines = format_otu_table(\
     sample_ids, otu_ids, otu_table_data, taxa).split('\n')
    
    result = filter_samples_from_otu_table(otu_table_lines,
                                           samples_to_keep,
                                           negate=True)
    return result
    
def filter_samples_from_distance_matrix(dm_lines,samples_to_discard,negate=False):
    """ Remove specified samples from distance matrix """
    sample_ids, dm_data = parse_distmat(dm_lines)
    
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
    