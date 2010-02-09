#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Justin Kuczynski","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"
 
from sys import argv, exit, stderr, stdout
from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.parse import fields_to_dict
from qiime.make_otu_table import make_otu_map

script_description = """
Makes sample x OTU table from OTU map and taxonomy.

Assumes that in the OTU map, the ids are in the format lib_seq, e.g.
M3FclSwb_1023. Will not work if this assumption is not met. Splits on last
underscore only so should be relatively robust to underscore in sample id."""

script_usage = """
Example 1:

Usage: make_otu_table.py -i seqs_otus.txt -o otu_table.txt

Example 2 - With a taxonomy file:

Usage: make_otu_table.py -i seqs_otus.txt -o otu_table.txt -t repr_set_tax_assignments.txt
"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
 make_option('-i', '--input_otu_fname', dest='otu_fname', help='Path to OTU \
file containing sequence ids assigned to each OTU (i.e., resulting OTU file \
from pick_otus.py)')
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
 make_option('-t', '--taxonomy', dest='taxonomy_fname', \
help='Path to taxonomy assignment, containing the assignments of taxons to \
sequences (i.e., resulting txt file from assign_taxonomy.py) \
[default: %default]', default=None),
 make_option('-o', '--output_fname', dest='output_fname', help='This is the \
filename that should be used when writing the output [default is stdout]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if opts.output_fname:
        outfile = open(opts.output_fname, 'w')
    else:
        outfile = stdout
    if not opts.taxonomy_fname:
        otu_to_taxonomy = None
    else:
        res = {}
        infile = open(opts.taxonomy_fname,'U')
        for line in infile:
            fields = line.split('\t')
            # typically this looks like: 3 SAM1_32 \t Root,Bacteria,Fi... \t 0.9
            # implying otu 3; sample 1, seq 32 (the representative of otu 3);
            # followed by the taxonomy and confidence
            if not len(fields) == 3:
                continue
            otu = fields[0].split(' ')[0]
            res[otu] = fields[1]
        otu_to_taxonomy = res

    otu_to_seqid = fields_to_dict(open(opts.otu_fname, 'U'))

    outfile.write(make_otu_map(otu_to_seqid, otu_to_taxonomy))
    

if __name__ == "__main__":
    main()