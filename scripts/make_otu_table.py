#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Justin Kuczynski","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "0.92"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"
 
from sys import argv, exit, stderr, stdout
from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from qiime.parse import fields_to_dict, parse_taxonomy
from qiime.make_otu_table import make_otu_map

options_lookup = get_options_lookup()

#make_otu_table.py
script_info={}
script_info['brief_description']="""Make OTU table"""
script_info['script_description']="""The script make_otu_table.py tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU in the last column if a taxonomy file is supplied."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""For this example the input is an OTU file containing sequence ids assigned to each OTU (i.e., resulting OTU file from pick_otus.py) and a text file containing the taxonomy assignments (i.e., resulting text file from assign_taxonomy.py), where the output file is defined as otu_table.txt:""","""%prog -i seqs_otus.txt -t repr_set_tax_assignments.txt -o otu_table.txt"""))
script_info['output_description']="""The output of make_otu_table.py is a tab-delimited text file, where the columns correspond to Samples and rows correspond to OTUs and the number of times a sample appears in a particular OTU."""
script_info['required_options']=[\
 options_lookup['otu_map_as_primary_input']
]
script_info['optional_options']=[ \
  make_option('-t', '--taxonomy', dest='taxonomy_fname', \
              help='Path to taxonomy assignment, containing the assignments of \ taxons to sequences (i.e., resulting txt file from assign_taxonomy.py) \
 [default: %default]', default=None),
  options_lookup['output_fp']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.output_fp:
        outfile = open(opts.output_fp, 'w')
    else:
        outfile = stdout
    if not opts.taxonomy_fname:
        otu_to_taxonomy = None
    else:
       infile = open(opts.taxonomy_fname,'U')
       otu_to_taxonomy = parse_taxonomy(infile)

    otu_to_seqid = fields_to_dict(open(opts.otu_map_fp, 'U'))

    outfile.write(make_otu_map(otu_to_seqid, otu_to_taxonomy))
    

if __name__ == "__main__":
    main()
