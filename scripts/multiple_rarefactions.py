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
import os.path
from qiime.rarefaction import RarefactionMaker

script_description = """Description:
subsample an otu table without replacement, to generate a collection of new
rarefied otu tables"""

script_usage = """
 Subsample otu_table.txt:
 python %prog otu_table.txt -m 100 -x 1200 -s 100 -n 2 -o mo/rare --small_included
(subsample otu_table.txt w/o replacement at 100 seqs per sample (twice),
200 seqs per sample (twice) ... 1200 seqs per sample (twice).
write 24 files total to mo/rare directory, named e.g.:rarefaction_100_0.txt
(100 seqs/sample; iteration 0, or first file written at 100 seqs/sam)
without --small_included, 24 output files are not guaranteed

 Subsample otu_table.txt repeatedly, at only one depth of seqs/sample:
 python %prog otu_table.txt -m 2500 -x 2500 -s 100 -n 100 -o mo/rare2 --small_included
100 files will be written to mo/rare2
 """

required_options = [
    make_option('-i', '--input_path',
        help='input otu table filepath'),

    make_option('-o', '--output_path',
        help='write output rarefied otu tables files to this dir ' +\
        "makes dir if it doesn't exist"),

    make_option('-m', '--min', type=int,
    help='min seqs/sample'),
    
    make_option('-x', '--max', type=int,
    help='max seqs/sample (inclusive)'),
    
    make_option('-s', '--step', type=int,
    help='levels: min, min+step... for level <= max')
]

optional_options = [

    make_option('-n', '--num-reps', dest='num_reps', default=1, type=int,
        help='num iterations at each seqs/sample level [default: %default]'),
     
    make_option('--lineages_included', dest='lineages_included', default=False,
        action="store_true",
          help="""output rarefied otu tables will include taxonomic (lineage)
          information for each otu, if present in input otu table
          [default: %default]"""),

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

    if opts.step <= 0:
        option_parser.error("nonpositive step not allowed (%s was supplied)" % \
          (opts.step,))
    if not os.path.exists(opts.output_path):
        os.makedirs(opts.output_path)
    maker = RarefactionMaker(opts.input_path, opts.min, opts.max,
        opts.step, opts.num_reps)
    maker.rarefy_to_files(opts.output_path, opts.small_included,
        include_lineages=opts.lineages_included)


if __name__ == "__main__":
    main()