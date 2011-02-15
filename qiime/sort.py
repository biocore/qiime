#!/usr/bin/env python
# File created on 15 Feb 2011

from __future__ import division
import re
from operator import itemgetter
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import (parse_mapping_file)

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso",
               "Rob Knight"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


def _natsort_key(item):
    """Provides normalized version of item for sorting with digits.

    From: 
    http://lists.canonical.org/pipermail/kragen-hacks/2005-October/000419.html
    """
    chunks = re.split('(\d+(?:\.\d+)?)', item)
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