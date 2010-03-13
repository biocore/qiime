#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Justin Kuczynski","Jeremy Widmann"] 
#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.92-dev"
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
    cols: metrics/ confidence levels, etc

    data: array or 2d list.
    """
    len_col = len(col_names)
    try:
        if data.shape != (len(row_names), len_col):
            raise ValueError, "Data shape of %s doesn't match header sizes %s %s" %\
                (data.shape, len(row_names), len(col_names))
    except AttributeError:
        # must be list of list
        try:
            if not numpy.all([len_col==len(row) for row in data]) or\
                    len(row_names) != len(data):
                raise ValueError, "Data shape doesn't match header sizes %s %s" %\
                    (len(row_names), len(col_names))
        except:
            raise ValueError, "Unsupported data type for format_matrix"

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
    """Writes OTU table as tab-delimited text.
    
    inputs: sample_names, otu_names are lists of strings
    data is numpy 2d array, num_otus x num_samples
    taxonomy is list of length = num_otus
    
    """
    if data.shape != (len(otu_names), len(sample_names)):
        raise ValueError, "Data shape of %s doesn't match %s OTUs, %s samples" \
            % (data.shape, len(otu_names), len(sample_names))
    lines = []
    #data = numpy.array(data, dtype='str') ##BAD! truncates some ints!
    sample_names = map(str, sample_names)
    otu_names = map(str, otu_names)
    lines.append('#'+comment)
    if taxonomy:
        lines.append('\t'.join(['#OTU ID'] + sample_names + 
            ['Consensus Lineage']))
        for otu_name, vals, taxon in zip(otu_names, data, taxonomy):
            if not isinstance(taxon, str):
                taxon = ';'.join(taxon)
            lines.append('\t'.join([otu_name] + map(str, vals.tolist()) + [taxon]))
    else:
        lines.append('\t'.join(['#OTU ID'] + sample_names))
        for otu_name, vals in zip(otu_names, data):
            lines.append('\t'.join([otu_name] + map(str,vals.tolist())))
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

def build_prefs_string(color_by_string):
    if not color_by_string:
        return ''
    fields = color_by_string.split(',')
    l = ['{']
    first = True
    entry_string = \
     "\t'%s':\n\t{\n\t\t'column':'%s',\n\t\t'colors':(('red',(0,100,100)),('blue',(240,100,100)))\n\t}"
    for field in fields:
        if first:
            first=False
            l.append('\n')
        else:
            l.append(',\n')
        l.append(entry_string % (field, field))
    l.append('\n}')
    return ''.join(l)

def format_map_file(headers, id_map, desc_key, sample_id_key, \
    description_map=None, run_description=None):
    """Generates string for formatted map file.
    
    Input:
        headers: list of strings corresponding to col headers
        id_map: dict of {id:{header:val}}
        description_map: dict of {id:description}
        run_description: either string, or list of strings
    """
    result = []
    if desc_key in headers:
        headers.remove(desc_key)
    if sample_id_key in headers:
        headers.remove(sample_id_key)
    header_line = '\t'.join([sample_id_key] + headers + [desc_key])
    if not header_line.startswith('#'):
        header_line = '#' + header_line
    result.append(header_line)
    if run_description:
        if not isinstance(run_description, str):
            run_description = '\n#'.join(run_description)
        if not run_description.startswith('#'):
            run_description = '#'+run_description
        result.append(run_description)
    for id_, fields in sorted(id_map.items()):
        curr_line = [id_]
        curr_line.extend([fields.get(h,'') for h in headers])
        curr_line.append(description_map.get(id_,''))
        result.append('\t'.join(map(str, curr_line)))
    return '\n'.join(result)
    
def format_histograms(pre_hist, post_hist, bin_edges):
    """Returns text-formatted histogram."""
    lines = []
    lines.append('Length\tBefore\tAfter')
    for edge, pre, post in zip(bin_edges, pre_hist, post_hist):
        lines.append('\t'.join(map(str, [edge, pre, post])))
    return '\n'.join(lines)
