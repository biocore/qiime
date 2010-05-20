#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.rarefaction import SingleRarefactionMaker

script_info={}
script_info['brief_description']="""Perform rarefaction on an otu table"""
script_info['script_description']="""To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, an OTU table.  This does not provide curves of diversity by number of sequences in a sample. Rather it creates a subsampled OTU table by random sampling (without replacement) of the input OTU table.  Samples that have fewer sequences then the requested rarefaction depth are omitted from the ouput otu tables.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""subsample otu_table.txt at 400 seqs/sample (-d), write results to a file (i.e. rarefaction_400_17.txt) ""","""single_rarefaction.py -i otu_table.txt -o rarefaction_400_17.txt -d 400"""))
script_info['script_usage'].append(('',"""(naming convention rarefaction_400_17.txt implies that the depth is 400 seqs/sam, iteration 17 at that depth (18th file written, due to iter 0))""",''))
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

make_option('--lineages_included', dest='lineages_included', default=False,
    action="store_true",
    help="""output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: %default]"""),

    make_option('-k', '--keep_empty_otus', default=False, action='store_true',
        help='otus (rows) of all zeros are usually omitted from the output otu table, with -k they will not be removed from the output file [default: %default]'),

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    maker = SingleRarefactionMaker(opts.input_path, opts.depth)
    maker.rarefy_to_file(opts.output_path, False,
        include_lineages=opts.lineages_included,
        empty_otus_removed=(not opts.keep_empty_otus))

if __name__ == "__main__":
    main()