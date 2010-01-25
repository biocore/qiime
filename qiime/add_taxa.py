#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

"""Contains code for adding taxa to OTU table that lacks them.
"""
from sys import stdout, stderr
from optparse import OptionParser
from string import strip
from qiime.parse import fields_to_dict

def fix_taxonomy_delimiters(taxonomy):
    """fixes delimiters in taxonomy (expect semicolons, but get commas)"""
    result = {}
    for k, vals in taxonomy.iteritems():
        v = vals[0]
        if ';' in v:
            result[k] = v
        else:
            result[k] = v.replace(',',';').replace('"','')
    return result

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options]'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-O','--otu_file',action='store',\
          type='string',dest='otu_fp',help='Path to read '+\
          'otu file [required]')
           
    parser.add_option('-T','--taxonomy_file',action='store',\
          type='string',dest='taxon_fp',help='Path to read '+\
          'taxonomy file [required]')
           
    parser.add_option('-o','--output_file',action='store',\
          type='string',dest='out_fp',help='Path to write '+\
          'output file')

    parser.add_option('-I','--id_map_file',action='store',\
          type='string',dest='id_map_fp',help='Path to read '+\
          'seq id to otu map file')

    opts,args = parser.parse_args()

    return opts, args
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    output_fname = opts.out_fp
    #Note: the OTU table is often too large to read into memory, hence
    #working on file directly.
    if output_fname:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout
    
    taxonomy = fields_to_dict(open(opts.taxon_fp, 'U'))
    taxonomy = fix_taxonomy_delimiters(taxonomy)

    if opts.id_map_fp:
        id_map = dict([map(strip, line.split('\t')) for line in 
            open(opts.id_map_fp, 'U')])
        new_taxonomy = dict([(id_map[k], v) for k, v in taxonomy.items()
            if k in id_map])
        assert new_taxonomy != taxonomy
        taxonomy = new_taxonomy

    for line in open(opts.otu_fp, 'U'):
        if line.startswith('#OTU ID'):
            outfile.write(line[:-1]+'\tConsensus Lineage\n')
        elif line.startswith('#'):
            outfile.write(line)
        else:
            id_, rest = line.split('\t', 1)
            t = taxonomy.get(id_, 'None')
            outfile.write(line[:-1]+'\t'+t+'\n')
