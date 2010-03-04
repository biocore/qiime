#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.rarefaction import SingleRarefactionMaker

script_info={}
script_info['brief_description']="""Perform rarefaction on a single otu table"""
script_info['script_description']="""To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, an OTU table.  This does not provide curves of diversity by number of sequences in a sample. Rather it creates a subsampled OTU table by random sampling (without replacement) of the input OTU table.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""subsample otu_table.txt at 400 seqs/sample (-d), write results to a file (i.e. rarefaction_400_17.txt) (samples which have fewer than 400 sequences will be included without subsampleing, exactly as they appear in otu_table.txt""","""single_rarefaction.py -i otu_table.txt -o rarefaction_400_17.txt -d 400 --small_included"""))
script_info['script_usage'].append(('',"""(naming convention implies that the depth is 200 seqs/sam, iteration 17 at that depth (18th file written, due to iter 0))""",''))
script_info['output_description']="""The results of single_rarefaction.py consist of a single subsampled OTU table. The file has the same otu table format as the input otu_table.txt. note: if the output file would be empty, no file is written"""


script_info['required_options']=[
    make_option('-i', '--input_path',
        help='input otu table filepath'),
    make_option('-o', '--output_path',
        help='write output rarefied otu tables to this filepath'),
    make_option('-d', '--depth', type=int,
        help='sequences per sample to subsample'),
]
script_info['optional_options']=[
make_option('--small_included', dest='small_included', default=False,
    action="store_true",
    help="""samples containing fewer seqs than the rarefaction level are included in the output but not rarefied [default: %default]"""),

make_option('--lineages_included', dest='lineages_included', default=False,
    action="store_true",
    help="""output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: %default]"""),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    maker = SingleRarefactionMaker(opts.input_path, opts.depth)
    maker.rarefy_to_file(opts.output_path, opts.small_included,
        include_lineages=opts.lineages_included)

if __name__ == "__main__":
    main()