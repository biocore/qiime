#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

import numpy
from numpy import isnan
from StringIO import StringIO

"""Contains formatters for the files we expect to encounter in 454 workflow.

A lot of this might migrate into cogent at some point.
"""

def format_distance_matrix(labels, data):
    """Writes distance matrix as tab-delimited text, uses format_matrix"""
    return format_matrix(data, labels, labels)

def format_matrix(data, row_names, col_names):
    """Writes matrix as tab-delimited text.

    format is rows: samples, 
    cols: metrics/ confidence levels, etc"""
    if data.shape != (len(row_names), len(col_names)):
        raise ValueError, "Data shape of %s doesn't match header sizes %s %s" %\
            (data.shape, len(row_names), len(col_names))
    lines = []
    row_names = map(str, row_names)   
    col_names = map(str, col_names)   
    #just in case they weren't strings initially
    lines.append('\t'.join([''] + col_names))
    for sam, vals in zip(row_names, data):
        lines.append('\t'.join([sam] + map(str, vals)))
    return '\n'.join(lines)

def format_otu_table(sample_names, otu_names, data, taxonomy=None,
    comment='Full OTU Counts'):
    """Writes OTU table as tab-delimited text."""
    if data.shape != (len(otu_names), len(sample_names)):
        raise ValueError, "Data shape of %s doesn't match %s OTUs, %s samples" \
            % (data.shape, len(otu_names), len(sample_names))
    lines = []
    data = numpy.array(data, dtype='str')
    sample_names = map(str, sample_names)
    otu_names = map(str, otu_names)
    lines.append('#'+comment)
    if taxonomy:
        lines.append('\t'.join(['#OTU ID'] + sample_names + 
            ['Consensus Lineage']))
        for otu_name, vals, taxon in zip(otu_names, data, taxonomy):
            if not isinstance(taxon, str):
                taxon = ';'.join(taxon)
            lines.append('\t'.join([otu_name] + vals.tolist() + [taxon]))
    else:
        lines.append('\t'.join(['#OTU ID'] + sample_names))
        for otu_name, vals in zip(otu_names, data):
            lines.append('\t'.join([otu_name] + vals.tolist()))
    return '\n'.join(lines)

def format_coords(coord_header, coords, eigvals, pct_var):
    """formats coords given specified coords matrix etc."""
    result = []
    result.append('pc vector number\t' +
        '\t'.join(map(str, range(1,len(coord_header)+1))))
    for name, row in zip(coord_header, coords):
        result.append('\t'.join([name] + map(str, row)))
    result.append('')
    result.append('')
    result.append('eigvals\t' + '\t'.join(map(str,eigvals)))
    result.append('% variation explained\t' +
        '\t'.join(map(str, pct_var)))
    return '\n'.join(result)

