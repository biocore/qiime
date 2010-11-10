#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 
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
    otu_lookup = {}.fromkeys([e.split()[0] for e in otus_to_discard])
    result = []
    
    if negate:
        def keep_seq(s):
            return s in otu_lookup
    else:
        def keep_seq(s):
            return s not in otu_lookup
    
    for line in otu_table_lines:
        line = line.strip()
        if line:
            if line.startswith('#'):
                # keep all comment lines
                result.append(line)
            elif keep_seq(line.split()[0]):
                result.append(line)
            else:
                # this line is filtered
                pass
            
    return result
    