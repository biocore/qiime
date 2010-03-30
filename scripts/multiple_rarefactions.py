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
import os.path
from qiime.rarefaction import RarefactionMaker

script_info={}
script_info['brief_description']="""Perform multiple subsamplings/rarefactions on an otu table"""
script_info['script_description']="""To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, OTU tables.  This does not provide curves of diversity by number of sequences in a sample. Rather it creates a series of subsampled OTU tables by random sampling (without replacement) of the input OTU table.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""An example of this script, where the user sets the minimum ("-m") and maximum ("-x") number of sequences per sample to 100 and 1200, respectively, while using steps ("-s") of 100, performing 2 iterations at each sampling depth ("-n"), and outputting the results to the directory "rarefaction_tables/" is shown by the following command:""","""multiple_rarefactions.py otu_table.txt -m 100 -x 1200 -s 100 -n 2 -o rarefaction_tables/"""))
script_info['script_usage'].append(("""""","""As a result, this command produces subsamples of the input otu_table.txt at 100 seqs per sample (twice), 200 seqs per sample (twice) ... 1200 seqs per sample (twice), which produces 24 rarefied otu talbes in the "rarefaction_tables" directory.""",""""""))
script_info['script_usage'].append(("""""","""Any sample containing fewer sequences in the input file than the requested number of sequences per sample is removed from the output rarefied otu table. To include samples with fewer than the requested number, you must manually add those samples to the resulting otu tables""",""""""))
script_info['output_description']="""The result of multiple_rarefactions.py consists of a number of files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.txt, and are named in the following way: rarefaction_100_0.txt, where "100" corresponds to the sequences per sample and "0" the iteration."""
script_info['required_options']=[
    make_option('-i', '--input_path',
        help='input otu table filepath'),

    make_option('-o', '--output_path',
        help="write output rarefied otu tables files to this dir makes dir if it doesn't exist"),

    make_option('-m', '--min', type=int,
    help='min seqs/sample'),
    
    make_option('-x', '--max', type=int,
    help='max seqs/sample (inclusive)'),
    
    make_option('-s', '--step', type=int,
    help='levels: min, min+step... for level <= max')
]


script_info['optional_options']=[

    make_option('-n', '--num-reps', dest='num_reps', default=1, type=int,
        help='num iterations at each seqs/sample level [default: %default]'),
     
    make_option('--lineages_included', dest='lineages_included', default=False,
        action="store_true",
          help="""output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: %default]"""),

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.step <= 0:
        option_parser.error("nonpositive step not allowed (%s was supplied)" % \
          (opts.step,))
    if not os.path.exists(opts.output_path):
        os.makedirs(opts.output_path)
    maker = RarefactionMaker(opts.input_path, opts.min, opts.max,
        opts.step, opts.num_reps)
    maker.rarefy_to_files(opts.output_path, False,
        include_lineages=opts.lineages_included)


if __name__ == "__main__":
    main()