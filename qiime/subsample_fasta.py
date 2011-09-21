#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@gmail.com"
__status__ = "Development"

from random import random

from cogent.parse.fasta import MinimalFastaParser


""" library functions for randomly subsampling fasta files.

"""

def subsample_fasta(input_fasta_fp,
                    output_fp,
                    percent_subsample):
    """ Writes random percent_sample of sequences from input fasta filepath
    
    input_fasta_fp: input fasta filepath
    output_fp: output fasta filepath
    percent_subsample: percent of sequences to write
    """
    
    input_fasta = open(input_fasta_fp, "U")
    
    output_fasta = open(output_fp, "w")

    for label, seq in MinimalFastaParser(input_fasta):
        if random() < percent_subsample:
            output_fasta.write('>%s\n%s\n' % (label, seq))
            

    
