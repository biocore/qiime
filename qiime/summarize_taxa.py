#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
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

    opts,args = parser.parse_args()

    return opts, args
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    output_fname = opts.out_fp
    result = {}
    #Note: the OTU table is often too large to read into memory, hence
    #working on file directly.
    if output_fname:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout
    level = opts.level 
    for line in open(opts.otu_fp, 'U'):
        if line.startswith('#OTU ID'):
            outfile.write(line.replace('#OTU ID', 'Taxon'))
        elif line.startswith('#'):
            outfile.write(line)
        else:
            fields = line.split('\t')
            data = array(fields[1:-1], dtype=int)
            tax = map(strip, fields[-1].split(';'))[:level]
            if len(tax) < level:
                tax.append('Other')
            tax_str = ';'.join(tax)
            if tax_str in result:
                result[tax_str] += data
            else:
                result[tax_str] = data
    for key, val in sorted(result.items()):
        outfile.write('\t'.join([key] + map(str, val))+'\n')
