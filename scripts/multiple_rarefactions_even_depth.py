#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.util import parse_command_line_parameters
from qiime.util import make_option
import os.path
from qiime.rarefaction import RarefactionMaker

script_info = {}
script_info[
    'brief_description'] = """Perform multiple rarefactions on a single otu table, at one depth of sequences/sample"""
script_info[
    'script_description'] = """To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, an OTU table.  This does not provide curves of diversity by number of sequences in a sample.  Rather it creates a subsampled OTU table by random sampling (without replacement) of the input OTU table.  Samples that have fewer sequences then the requested rarefaction depth are omitted from the ouput otu tables.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG."""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Example:""",
     """subsample otu_table.biom at 100 seqs/sample (-d) 10 times (-n) and write results to files (e.g., rarefaction_400_0.biom) in 'rarefied_otu_tables/' (-o).""",
     """%prog -i otu_table.biom -o rarefied_otu_tables/ -d 100 -n 10"""))

script_info[
    'output_description'] = """The results of this script consist of n subsampled OTU tables, written to the directory specified by -o. The file has the same otu table format as the input otu_table.biom. Note: if the output files would be empty, no files are written."""


script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_filepath',
                help='input otu table filepath'),

    make_option('-o', '--output_path', type='new_dirpath',
                help="write output rarefied otu tables files to this dir (makes dir if it doesn't exist)"),

    make_option('-d', '--depth', type='int',
                help='sequences per sample to subsample'),
]


script_info['optional_options'] = [

    make_option('-n', '--num_reps', dest='num_reps', default=10, type='int',
                help='num iterations at each seqs/sample level [default: %default]'),

    make_option('--lineages_included', dest='lineages_included', default=False,
                action="store_true",
                help="""output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: %default]"""),

    make_option('-k', '--keep_empty_otus', default=False, action='store_true',
                help='otus (rows) of all zeros are usually omitted from the output otu tables, with -k they will not be removed from the output files [default: %default]'),

]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if not os.path.exists(opts.output_path):
        os.makedirs(opts.output_path)
    maker = RarefactionMaker(opts.input_path, opts.depth, opts.depth,
                             1, opts.num_reps)
    maker.rarefy_to_files(opts.output_path, False,
                          include_lineages=opts.lineages_included,
                          empty_otus_removed=(not opts.keep_empty_otus))


if __name__ == "__main__":
    main()
