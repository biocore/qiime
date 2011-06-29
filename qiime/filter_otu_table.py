#!/usr/bin/env python
#filter_otu_table
"""Filters OTU table according to minimum OTU count and number of samples.

If OTU table has taxonomy assigned, can also use taxonomy to filter.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Rob Knight and Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Tony Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Release"

from qiime.parse import parse_otu_table
from qiime.format import format_otu_table
from qiime.util import create_dir
from string import strip
from numpy import array
import os


def strip_quotes(s):
    """splits leading/trailing quotes from string s"""
    if not s or len(s) < 2:
        return s
    if s[0] == s[-1]:
        if s[0] in '"\'':
            s = s[1:-1]
    return s

def split_tax(tax):
    """splits tax string on semicolon and comma"""
    fields = tax.split(';')
    if len(fields) == 1:
        fields = fields[0].split(',')
    return map(strip_quotes, fields)

def filter_table(params,filtered_table_path,otu_file):
    """ Filters table according to OTU counts, occurance, and taxonomy
    
    params: Dictionary containing minimum sequence count (min_otu_count) 
     per OTU, minimum number of samples that OTU needs to occur in
     (min_otu_samples), targetted taxonomy to retain (included_taxa), and
     taxonomy to exclude (excluded_taxa). 
    filtered_table_path:  Open file object to write filtered table to.
    otu_file: Open file object of input OTU file.
    """
    

    min_otu_count=params['min_otu_count']
    min_otu_samples=params['min_otu_samples']
    included_taxa=params['included_taxa']
    excluded_taxa=params['excluded_taxa']
    
    otu_data = parse_otu_table(otu_file)
    
    # Create list of OTUs that fail to pass filters
    flagged_otus = []
    
    otu_index = 1
    otus = otu_data[otu_index]
    
    otu_counts_index = 2
    otu_counts = otu_data[otu_counts_index]
    
    taxa_index = 3
    
    try:
        taxa_lines = otu_data[taxa_index]
        if len(taxa_lines):
            taxa_present = True
        else:
            taxa_present = False
    except IndexError:
        taxa_present = False
    
    index_counter = -1
    for otu_count in otu_counts:

        index_counter += 1
        
        if otu_count.sum() < min_otu_count or \
         (otu_count > 0).sum() < min_otu_samples:
             flagged_otus.append(otus[index_counter])
             continue
        if taxa_present:
            taxa = set(taxa_lines[index_counter])
            # Check for targetted taxa that also are not excluded
            if taxa.intersection(included_taxa) and not \
             taxa.intersection(excluded_taxa):
                continue
            # If taxonomy found in included taxa and no excluded taxa 
            # given, skip filtering.
            elif taxa.intersection(included_taxa) and not excluded_taxa:
                continue
            # Skip any taxonomic filtering if taxa present but no filters given
            elif not included_taxa and not excluded_taxa:
                continue
            # If only specifying exluded taxa, allow inclusion of this OTU
            # if taxa not in excluded set.
            elif not included_taxa and not taxa.intersection(excluded_taxa):
                continue
            # taxa does is not included, or falls in excluded set, so flag
            # this OTU for removal
            else:
                flagged_otus.append(otus[index_counter])

    sample_id_index = 0
        
    raw_otu_table = (format_otu_table(otu_data[sample_id_index], 
     otus, otu_counts, taxonomy=taxa_lines, skip_empty=True)).split('\n')
     
    # Filter out lines of the OTU table that are flagged
    
    filtered_otu_table = ""
    
    for line in raw_otu_table:
        if line.startswith("#"):
            filtered_otu_table += line + '\n'
            continue
        curr_otu_id = line.split('\t')[0].strip()
        
        if curr_otu_id in flagged_otus:
            continue
        else:
            filtered_otu_table += line + '\n'
        
    filtered_table_path.write(filtered_otu_table)
  
    

def _filter_table_samples(otu_table_lines, min_seqs_per_sample):
    """removes samples from OTU_table that have less than min_seqs_per_sample
    """
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(otu_table_lines)
    counts = sum(otu_table)
    big_enough_samples = (counts>=int(min_seqs_per_sample)).nonzero()
    res_otu_table = otu_table.copy()
    res_otu_table = res_otu_table[:,big_enough_samples[0]]
    res_sample_ids = map(sample_ids.__getitem__, big_enough_samples[0])
    return format_otu_table(res_sample_ids, otu_ids, res_otu_table, lineages)


def _filter_table_neg_control(otu_table_lines, samples):
    """removes OTUs from OTU_table that are found in one of the samples in the sample list
    """
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(otu_table_lines)
    new_otu_table = []
    new_otu_ids = []
    new_lineages = []
    #get the sample indices to remove
    sample_indices = []
    for i in samples:
        if i in sample_ids:
            index = sample_ids.index(i)
            sample_indices.append(index)

    for i, row in enumerate(otu_table):
        #figure out if the OTU is in any of the negative controls
        count = 0
        for j in sample_indices:
            count += row[j]
        #only write it to the new OTU table if it is not
        if count == 0:
            if lineages:
                new_lineages.append(lineages[i])
            new_otu_table.append(list(row))
            new_otu_ids.append(otu_ids[i])
    new_otu_table = array(new_otu_table)
    result = format_otu_table(sample_ids, new_otu_ids, new_otu_table, new_lineages)
    result = result.split('\n')
    #remove the samples
    return _filter_table_samples(result, 1)


