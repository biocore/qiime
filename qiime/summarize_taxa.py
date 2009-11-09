#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight", "Catherine Lozupone"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

"""Contains code for summarizing OTU table with taxa in last field.
"""
from sys import stdout, stderr
from optparse import OptionParser
from string import strip
from numpy import array
from otu_category_significance import convert_OTU_table_relative_abundance

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options]'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-O','--otu_file',action='store',\
          type='string',dest='otu_fp',help='Path to read '+\
          'otu file [required]')
           
    parser.add_option('-o','--output_file',action='store',\
          type='string',dest='out_fp',help='Path to write '+\
          'output file')
    
    parser.add_option('-L','--level',action='store',\
          type='int',dest='level', default=2, 
          help='Level of taxonomy to use')

    parser.add_option('-m','--category_mapping',action='store',\
          type='string',dest='category_mapping', 
          help='if supplied - the taxon information will be added ' +\
             'to the category mapping file. This mapping file can ' +\
             'be used to color PCoA plots by taxon abundance or to ' +\
             'perform statistical tests of taxon/category associations.')
   
    parser.add_option('-d','--delimitor',action='store',\
          type='string',dest='delimitor',default=';', 
          help='Delimitor that separates taxonomy categories. default=;')
    
    parser.add_option('-r', '--relative_abundance', action='store',\
        type='string', dest='relative_abundance', default='True', \
        help='If True, reports the relative abundance of the lineage in ' +\
            'each sample. If False, reports the raw counts. Default=True')

    opts,args = parser.parse_args()
    return opts, args

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
            output.append(line + '\n')
        elif line.startswith('#'):
            output.append(line + '\n')
        else:
            result = process_data_line(line, result, delimitor, level)
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
            result = process_data_line(line, result, delimitor, level)
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
            index = samples.index(sample_name)
            for key, val in sorted(result.items()):
                line.append(str(val[index]))
            output.append('\t'.join(line) + '\n')
    return output

def process_data_line(line, result, delimitor, level):
    """
    """
    fields = line.split('\t')
    vals = [float(i) for i in fields[1:-1]]
    data = array(vals, dtype=float)
    tax = map(strip, fields[-1].split(delimitor)[:level])
    if len(tax) < level:
        tax.append('Other')
    tax_str = ';'.join(tax)
    if tax_str in result:
        result[tax_str] += data
    else:
        result[tax_str] = data
    return result

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    output_fname = opts.out_fp
    otu_fp = opts.otu_fp
    otu_table = open(otu_fp, 'U')
    delimitor = opts.delimitor
    category_mapping = opts.category_mapping
    if category_mapping:
        category_mapping = open(category_mapping)
    level = opts.level 
    relative_abundance=opts.relative_abundance
    if output_fname:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout
    if not category_mapping:
        output = make_new_summary_file(otu_table, level, delimitor, \
            relative_abundance) 
    else:
        output = add_summary_category_mapping(otu_table, category_mapping, \
            level, delimitor, relative_abundance) 
    outfile.write(''.join(output))

