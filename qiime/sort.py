#!/usr/bin/env python
# File created on 15 Feb 2011

from __future__ import division
import re
from operator import itemgetter
from numpy import array
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import (parse_mapping_file, parse_otu_table)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso",
               "Rob Knight"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


def _natsort_key(item):
    """Provides normalized version of item for sorting with digits.

    From: 
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    try:
        chunks = re.split('(\d+(?:\.\d+)?)', item)
    except TypeError:
        # if item is a tuple or list (i.e., indexable, but not a string)
        # work with the first element
        chunks = re.split('(\d+(?:\.\d+)?)', item[0])
    for ii in range(len(chunks)):
        if chunks[ii] and chunks[ii][0] in '0123456789':
            if '.' in chunks[ii]: numtype = float
            else: numtype = int
            # wrap in tuple with '0' to explicitly specify numbers come first
            chunks[ii] = (0, numtype(chunks[ii]))
        else:
            chunks[ii] = (1, chunks[ii])
    return (chunks, item)

def natsort(seq):
    """Sort a sequence of text strings in a reasonable order.

    From: 
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    alist = list(seq)
    alist.sort(key=_natsort_key)
    return alist


def sort_sample_ids_by_mapping_value(mapping_file,field,field_type_f=float):
    """ Return list of sample ids sorted by ascending value from mapping file
    """
    data, headers, comments = parse_mapping_file(mapping_file)

    try:
        column = headers.index(field)
    except ValueError:
        raise ValueError,\
         "Column (%s) not found in mapping file headers:\n %s" %\
         (field,' '.join(headers))

    results = [(e[0],field_type_f(e[column])) for e in data]
    results.sort(key=itemgetter(1))
    return results
    
def sort_fasta_by_abundance(fasta_lines,fasta_out_f):
    """ Sort seqs in fasta_line by abundance, write all seqs to fasta_out_f

     Note that all sequences are written out, not just unique ones.

     fasta_lines: input file handle (or similar object)
     fasta_out_f: output file handle (or similar object)

    ** The current implementation works well for fairly large data sets, 
       (e.g., several combined 454 runs) but we may want to revisit if it
       chokes on very large (e.g., Illumina) files. --Greg **

    """
    seq_index = {}
    count = 0
    for seq_id,seq in MinimalFastaParser(fasta_lines):
        count += 1
        try:
            seq_index[seq].append(seq_id)
        except KeyError:
            seq_index[seq] = [seq_id]

    seqs = []
    for k,v in seq_index.items():
        seqs.append((len(v),k,v))
        del seq_index[k]
    seqs.sort()
    for count, seq, seq_ids in seqs[::-1]:
        for seq_id in seq_ids:
            fasta_out_f.write('>%s\n%s\n' % (seq_id,seq))

def sort_otu_table_by_mapping_field(otu_table_data,
                                    mapping_file_data,
                                    sort_field,
                                    sort_f=natsort):
    """ sort otu table based on the value of sort_field for each sample 
    """
    mapping_data, header_data, comments = mapping_file_data

    mapping_field_index = header_data.index(sort_field)
    sorted_sample_ids = [(e[mapping_field_index], e[0]) for e in mapping_data]
    sorted_sample_ids = sort_f(sorted_sample_ids)
    sorted_sample_ids = [e[1] for e in sorted_sample_ids]
    
    return sort_otu_table(otu_table_data,sorted_sample_ids)
    
def sort_otu_table(otu_table_data,
                   sorted_sample_ids):
    sample_ids, otu_ids, otu_data, taxa = otu_table_data
    unsorted_sample_id_index = dict([(e,i) for i,e in enumerate(sample_ids)])
    
    sorted_sample_ids_set = set(sorted_sample_ids)
    if set(sample_ids) - sorted_sample_ids_set:
        raise KeyError,\
         "Sample IDs present in OTU table but not sorted sample id list: " + \
         ' '.join(list(set(sample_ids) - set(sorted_sample_ids)))
    if len(sorted_sample_ids_set) != len(sorted_sample_ids):
        raise ValueError,\
         "Duplicate sample IDs are present in sorted sample id list."

    otu_data = otu_data.transpose()
    sample_id_results = []
    sorted_otu_data = []
    for current_sid in sorted_sample_ids:
        try:
            current_otu_data = otu_data[unsorted_sample_id_index[current_sid]]
        except KeyError:
            # ignore samples in mapping that are not in otu table
            continue
        sample_id_results.append(current_sid)
        sorted_otu_data.append(current_otu_data)
    sorted_otu_data = array(sorted_otu_data)
    sorted_otu_data = sorted_otu_data.transpose()

    return sample_id_results, otu_ids, sorted_otu_data, taxa
