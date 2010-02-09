#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from sys import stdout
from qiime.add_taxa import rewrite_otu_table_with_taxonomy

script_description = """Adds taxa to OTU table that lacks them. """

script_usage = """Add taxa to otu file from otus.txt from file taxa.txt, writing result to stdout:

add_taxa.py -i otus.txt -t taxa.txt
"""

required_options = [\
    make_option('-i','--otu_file',action='store',\
        type='string',dest='otu_fp',help='Path to read otu file'),
    make_option('-t','--taxonomy_file',action='store',\
        type='string',dest='taxon_fp',help='Path to read taxonomy file')
]

optional_options = [\
    make_option('-o','--output_file',action='store',\
        type='string',dest='out_fp',help='Path to write '+\
        'output file [default: stdout]'),
    make_option('-m','--id_map_file',action='store',\
        type='string',dest='id_map_fp',help='Path to read '+\
        'seq id to otu map file [default: %default]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if opts.out_fp:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout
    taxon_lines = open(opts.taxon_fp, 'U')
    if opts.id_map_fp:
        id_map_lines = open(opts.id_map_fp, 'U')
    else:
        id_map_lines = None

    otu_lines = open(opts.otu_fp, 'U')
    rewrite_otu_table_with_taxonomy(taxon_lines, otu_lines, id_map_lines,
        outfile)

if __name__ == "__main__":
    main()
