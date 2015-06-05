#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger",
               "Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"

from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.pick_rep_set import (rep_set_picking_methods,
                                reference_rep_set_picking_methods)

script_info = {}
script_info['brief_description'] = """Pick representative set of sequences"""
script_info['script_description'] = """After picking OTUs, you can then pick a\
 representative set of sequences. For each OTU, you will end up with one\
 sequence that can be used in subsequent analyses."""
script_info['script_usage'] = []

script_info['script_usage'].append(("""Simple example: picking a representative\
 set for de novo-picked OTUs""", """The script pick_rep_set.py takes as input\
 an 'OTU map' (via the \"-i\" parameter) which maps OTU identifiers to sequence\
 identifiers. Typically, this will be the output file provided by pick_otus.py.\
 Additionally, a FASTA file is required, via \"-f\", which contains all of the\
 sequences whose identifiers are listed in the OTU map.""",
     """%prog -i seqs_otus.txt -f seqs.fna -o rep_set1.fna"""))

script_info['script_usage'].append(
    ("""Picking OTUs with "preferred representative" sequences""",
     """Under some circumstances you may have a fasta file of "preferred\
 representative" sequences. An example of this is if you were to pick OTUs\
 against a reference collection with uclust_ref. In this case you may want your\
 representative sequences to be the sequences from the reference collection,\
 rather than the sequences from your sequencing run. To achieve this, you can\
 pass the original reference collection via -r. If you additionally allowed for\
 new clusters (i.e., sequences which don't match a reference sequence are used\
 as seeds for new OTUs) you'll also need to pass the original sequence\
 collection to pick a representative sequence from the sequencing run in that\
 case.""",
     """%prog -i seqs_otus.txt -f seqs.fna -r refseqs.fasta -o rep_set2.fna"""))

script_info['output_description'] = """The output from pick_rep_set.py is a\
 single FASTA file containing one sequence per OTU. The FASTA header lines will\
 be the OTU identifier (from here on used as the unique sequence identifier)\
 followed by a space, followed by the sequence identifier originally associated\
 with the representative sequence. The name of the output FASTA file will be\
 <input_sequences_filepath>_rep_set.fasta by default, or can be specified via\
 the \"-o\" parameter.
"""

script_info['required_options'] = [
    make_option('-i', '--input_file', action='store',
                type='existing_filepath', dest='otu_fp', help='Path to input otu '
                'mapping file [REQUIRED]')
]
rep_set_picking_method_choices = rep_set_picking_methods.keys()
script_info['optional_options'] = [
    make_option('-f', '--fasta_file', action='store',
                type='existing_filepath', dest='fasta_fp', help='Path to input '
                'fasta file [REQUIRED if not picking against a '
                'reference set; default: None]'),
    make_option('-m', '--rep_set_picking_method',
                type='choice', dest='rep_set_picking_method',
                help=('Method for picking representative sets.  Valid choices '
                      'are ' + ', '.join(rep_set_picking_method_choices) +
                      ' [default: %default (first chooses cluster seed when picking '
                      'otus with uclust)]'),
                choices=rep_set_picking_method_choices, default='first'),
    make_option('-o', '--result_fp', action='store',
                type='new_filepath', dest='result_fp', help='Path to store '
                'result file [default: <input_sequences_filepath>_rep_set.fasta]'),
    make_option('-l', '--log_fp', action='store',
                type='new_filepath', dest='log_fp', help='Path to store '
                'log file [default: No log file created.]'),
    make_option('-s', '--sort_by', action='store',
                type='choice', choices=['otu', 'seq_id'],
                dest='sort_by', default='otu',
                help='sort by otu or seq_id [default: %default]'),
    make_option('-r', '--reference_seqs_fp', type='existing_filepath',
                help='collection of preferred representative '
                'sequences [default: %default]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    reference_seqs_filepath = opts.reference_seqs_fp
    input_seqs_filepath = opts.fasta_fp
    input_otu_filepath = opts.otu_fp
    result_path = opts.result_fp or\
        '%s_rep_set.fasta' % input_seqs_filepath
    log_path = opts.log_fp

    if reference_seqs_filepath:
        rep_set_picker =\
            reference_rep_set_picking_methods[opts.rep_set_picking_method]
        rep_set_picker(input_seqs_filepath,
                       input_otu_filepath,
                       reference_seqs_filepath,
                       result_path=result_path,
                       log_path=log_path,
                       sort_by=opts.sort_by)
    else:
        if not input_seqs_filepath:
            option_parser.error('--fasta_fp must be provided when not picking'
                                ' representative against a reference set.')
        rep_set_picker =\
            rep_set_picking_methods[opts.rep_set_picking_method]
        rep_set_picker(input_seqs_filepath,
                       input_otu_filepath,
                       result_path=result_path,
                       log_path=log_path,
                       sort_by=opts.sort_by)


if __name__ == "__main__":
    main()
