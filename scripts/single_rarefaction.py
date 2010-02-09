#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.rarefaction import SingleRarefactionMaker

script_description = """Description:
subsample an otu table without replacement, to generate a new rarefied otu table"""

script_usage = """ 
 Subsample otu_table.txt at 20 seqs/sample.  if any samples
have fewer than 20 sequences, include them as they appear in otu_table.txt 
(don't subsample them).
  python %prog otu_table.txt -o rarefaction_20_17.txt -d 20 --small_included
(naming convention implies that the depth is 20 seqs/sam, iteration 17 at that
  depth (18th file written, due to iter 0))
  
note:
if the output file would be empty, no file is written"""
required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
 make_option('-i', '--input_path',
     help='input otu table filepath'),
 make_option('-o', '--output_path',
     help='write output rarefied otu tables to this filepath) '),
 make_option('-d', '--depth', type=int,
         help='sequences per sample to subsample'),
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
     make_option('--small_included', dest='small_included', default=False,
         action="store_true",
         help="""samples containing fewer seqs than the rarefaction
 level are included in the output but not rarefied [default: %default]""")
]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    maker = SingleRarefactionMaker(opts.input_path, opts.depth)
    maker.rarefy_to_file(opts.output_path, opts.small_included)

if __name__ == "__main__":
    main()