#!/usr/bin/env python
# File created on 08 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

def ascii_to_phred(c,offset):
    """  Convert ascii character to Phred quality score with specified ASCII offset.
    """
    return ord(c) - offset
    
def ascii_to_phred33(c):
    """  Convert ascii character to Phred quality score with ASCII offset of 33.
      Standard "Sanger" ASCII offset of 33. This is used by Illumina in CASAVA versions 
       after 1.8.0, and most other places. Note that internal Illumina files still
       use offset of 64
    """
    return ascii_to_phred(c,33)

def ascii_to_phred64(c):
    """ Convert ascii character to Phred quality score with ASCII offset of 64.
      Illumina-specific ASCII offset of 64. This is used by Illumina in CASAVA versions 
       prior to 1.8.0, and in Illumina internal formats (e.g., export.txt files).
    """
    return ascii_to_phred(c,64)

def phred_to_ascii(d,offset):
    """ Convert Phred quality score to ASCII character with specified offset.
    """
    return chr(d + offset)

def phred_to_ascii33(d):
    """ Convert Phred quality score to ASCII character with offset of 33.
    """
    return phred_to_ascii(d,33)
    
def phred_to_ascii64(d):
    """ Convert Phred quality score to ASCII character with offset of 64.
    """
    return phred_to_ascii(d,64)


