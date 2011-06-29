#!/usr/bin/env python
# File created on 27 Apr 2010
from __future__ import division

""" Sort denoiser output by cluster size"""

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Release"
 
from qiime.util import make_option
from qiime.util import parse_command_line_parameters, get_options_lookup
from cogent.parse.fasta import MinimalFastaParser
from qiime.format import write_Fasta_from_name_seq_pairs
from qiime.denoise_wrapper import extract_cluster_size

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Sort denoiser output by cluster size."
script_info['script_description'] = "This scripts is used prior to OTU picking when combining several separately denoised data sets. It sorts the FASTA file by cluster size, such that the OTU pciker now which are the most likely the best OTU centroids."
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example Usage:""","""""","""sort_denoiser_output.py -f denoised_seqs.fasta -o denoised_seqs_sorted.fasta"""))

script_info['output_description']= "A standard FASTA file"
script_info['required_options'] = [\
    # Example required option
    options_lookup['input_fasta'],
    make_option('-o','--output_file',help='the output filename'),\
        ]

script_info['optional_options'] = [\
    # Example optional option
]

script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    seqs_w_cluster_size = [(extract_cluster_size(name), name, seq) for name,seq
                           in MinimalFastaParser(open(opts.input_fasta_fp))]

    seqs_w_cluster_size.sort(reverse=True)
    name_seqs = [(name,seq) for (cs,name,seq) in seqs_w_cluster_size]
    try:
        out_fh = open(opts.output_file,"w")
    except OSError:
        #re-raise slightly more informative
        raise OSError,"Could not write to file %s. Check permissions and file name" %\
            opts.output_file
    
    write_Fasta_from_name_seq_pairs(name_seqs, out_fh)

if __name__ == "__main__":
    main()
