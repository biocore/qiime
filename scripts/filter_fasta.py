#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from optparse import make_option
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import fields_to_dict
from qiime.filter import filter_fasta

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
 options_lookup['input_fasta'],
 make_option('-o','--output_fasta_fp',help='the output fasta filepath'),\
 make_option('-m','--seqs_to_keep_map',help='an OTU map where sequences ids are'
  'those which should be retained.'),\
]
script_info['optional_options'] = [\
 make_option('-n','--negate',help='keep only seqs that are not in seqs_to_keep [default: %default]',default=False,action='store_true')
]
script_info['version'] = __version__

def filter_fasta_fp(input_seqs_fp,output_seqs_fp,seqs_to_keep,negate=False):
    """Filter a fasta file to include only sequences listed in seqs_to_keep """
    
    input_seqs = MinimalFastaParser(open(input_seqs_fp,'U'))
    output_f = open(output_seqs_fp,'w')
    return filter_fasta(input_seqs,output_f,seqs_to_keep,negate)


def get_seqs_to_keep_lookup(seqs_to_keep_f):
    otu_map = fields_to_dict(seqs_to_keep_f)
    seqs_to_keep = []
    for seq_ids in otu_map.values():
        seqs_to_keep += seq_ids
    return {}.fromkeys(seqs_to_keep)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    seqs_to_keep_lookup = get_seqs_to_keep_lookup(open(opts.seqs_to_keep_map,'U'))
    filter_fasta_fp(opts.input_fasta_fp,
                    opts.output_fasta_fp,
                    seqs_to_keep_lookup,
                    opts.negate)


if __name__ == "__main__":
    main()