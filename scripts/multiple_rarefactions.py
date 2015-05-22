#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

import os.path
from functools import partial

from skbio.stats import subsample

from qiime.util import parse_command_line_parameters, create_dir
from qiime.util import make_option
from qiime.rarefaction import RarefactionMaker

script_info = {}
script_info[
    'brief_description'] = """Perform multiple subsamplings/rarefactions on an otu table"""
script_info[
    'script_description'] = """To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, OTU tables.  This does not provide curves of diversity by number of sequences in a sample.  Rather it creates a series of subsampled OTU tables by random sampling (without replacement) of the input OTU table.  Samples that have fewer sequences then the requested rarefaction depth for a given output otu table are omitted from those ouput otu tables.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG."""
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Generate rarefied OTU tables:""",
     """Generate rarefied OTU tables beginning with 10 (-m) sequences/sample through 140 (-x) sequences per sample in steps of of 10 (-s), performing 2 iterations at each sampling depth (-n). All resulting OTU tables will be written to 'rarefied_otu_tables' (-o). Any sample containing fewer sequences in the input file than the requested number of sequences per sample is removed from the output rarefied otu table.""",
     """%prog -i otu_table.biom -m 10 -x 140 -s 10 -n 2 -o rarefied_otu_tables/"""))

script_info[
    'output_description'] = """The result of multiple_rarefactions.py consists of a number of biom files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.biom, and are named in the following way: rarefaction_100_0.biom, where "100" corresponds to the sequences per sample and "0" the iteration."""

script_info['required_options'] = [
    make_option('-i', '--input_path',
                help='Input OTU table filepath.',
                type='existing_filepath'),
    make_option('-o', '--output_path',
                help="Output directory.",
                type='new_dirpath'),
    make_option('-m', '--min', type='int',
                help='Minimum number of seqs/sample for rarefaction.'),
    make_option('-x', '--max', type='int',
                help='Maximum number of seqs/sample (inclusive) for rarefaction. '),
    make_option('-s', '--step', type='int',
                help='Size of each steps between the min/max of' +
                ' seqs/sample (e.g. min, min+step... for level <= max).')
]
script_info['optional_options'] = [
    # This option is screwed up: dest should equal the long form parameter name,
    # but I'm not sure if we can do anything about it since '-' is not allowed
    # in variable names... Hmmm... Changing the long-form parameter name
    # would cause older parameter files not to work.
    make_option('-n', '--num_reps', dest='num_reps', default=10, type='int',
                help='The number of iterations at each step. [default: %default]'),
    make_option('--lineages_included', default=False,
                action="store_true",
                help='Retain taxonomic (lineage) information for each OTU. Note:' +
                ' this will only work if lineage information is in the input OTU' +
                ' table. [default: %default]'),
    make_option('-k', '--keep_empty_otus', default=False, action='store_true',
                help='Retain OTUs of all zeros, which are usually omitted from' +
                ' the output OTU tables. [default: %default]'),
    make_option('--subsample_multinomial', default=False, action='store_true',
                help='subsample using subsampling with replacement [default: %default]')
]

script_info['option_label'] = {'input_path': 'OTU table filepath',
                               'output_path': 'Output directory',
                               'min': 'Min # of seqs/sample',
                               'max': 'Max # of seqs/sample',
                               'step': 'Step size',
                               'num_reps': '# of iterations',
                               'lineages_included': 'Include lineages',
                               'keep_empty_otus': 'Retain empty OTUs'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.step <= 0:
        option_parser.error("step must be greater than 0")

    create_dir(opts.output_path, fail_on_exist=False)
    maker = RarefactionMaker(opts.input_path, opts.min, opts.max,
                             opts.step, opts.num_reps)

    if opts.subsample_multinomial:
        subsample_f = partial(subsample, replace=True)
    else:
        subsample_f = subsample

    maker.rarefy_to_files(opts.output_path,
                          False,
                          include_lineages=opts.lineages_included,
                          empty_otus_removed=(not opts.keep_empty_otus),
                          subsample_f=subsample_f)


if __name__ == "__main__":
    main()
