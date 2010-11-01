#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Development"

"""Contains code for summarizing OTU table with taxa in last field.
"""
from sys import stdout, stderr
from optparse import OptionParser
from string import strip
from numpy import array
from qiime.util import convert_OTU_table_relative_abundance

def make_new_summary_file(otu_table, level, delimitor, relative_abundance): 
    """makes a summary file with taxa on rows and samples in columns
    """
    #Note: the OTU table is often too large to read into memory, hence
    #working on file directly.
    result = {}
    output = []
    if relative_abundance == 'True':
        otu_table = convert_OTU_table_relative_abundance(otu_table)
    for line in otu_table:
        if line.startswith('#OTU ID'):
            line = line.replace('#OTU ID', 'Taxon')
            line = line.replace('\tConsensus Lineage', '')
            output.append(line.strip('\n') + '\n')
        elif line.startswith('#'):
            output.append(line.strip('\n') + '\n')
        else:
            result = process_data_line(line, result, delimitor, level,
                relative_abundance)
    for key, val in sorted(result.items()):
        output.append('\t'.join([key] + map(str, val))+'\n')
    return output

def add_summary_category_mapping(otu_table, category_mapping, \
    level, delimitor, relative_abundance): 
    """makes a summary file with taxa on rows and samples in columns
    """
    #Note: the OTU table is often too large to read into memory, hence
    #working on file directly.
    result = {}
    output = []
    if relative_abundance == 'True':
        otu_table = convert_OTU_table_relative_abundance(otu_table)
    for line in otu_table:
        if line.startswith('#OTU ID'):
            samples = line.strip().split('\t')[1:-1]
        if not line.startswith('#'):
            result = process_data_line(line, result, delimitor, level,
                relative_abundance)
    for line in category_mapping:
        if line.startswith('#SampleID'):
            header = line.strip().split('\t')
            header.extend(sorted(result.keys()))
            output.append('\t'.join(header) + '\n')
        elif line.startswith('#'):
            output.append(line + '\n')
        else:
            line = line.strip().split('\t')
            sample_name = line[0]
            try:
                index = samples.index(sample_name)
            except ValueError:
                pass # don't include this sample, as was not found in otu table
            else:
                for key, val in sorted(result.items()):
                    line.append(str(val[index]))
                output.append('\t'.join(line) + '\n')
    return output

def process_data_line(line, result, delimitor, level, use_float='True'):
    """ beware use_float is a STRING!!!
    """
    fields = line.split('\t')
    vals = [float(i) for i in fields[1:-1]]
    if use_float == 'True':
        data = array(vals, dtype=float)
    else:
        data = array(vals, dtype=int)
    tax = map(strip, fields[-1].split(delimitor)[:level])
    if len(tax) < level:
        tax.append('Other')
    tax_str = ';'.join(tax)
    if tax_str in result:
        result[tax_str] += data
    else:
        result[tax_str] = data
    return result

