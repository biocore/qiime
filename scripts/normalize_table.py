#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from qiime.util import parse_command_line_parameters, make_option
from qiime.normalize_table import normalize_CSS, normalize_DESeq, multiple_file_normalize_CSS, multiple_file_normalize_DESeq, algorithm_list

import rpy2.robjects as robjects
import os

__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"



script_info = {}
script_info['brief_description'] = """Alternate (not rarefying) matrix normalization techniques"""
script_info['script_description'] = """To perform any downstream analysis after OTU picking (besides correlations and CSS/DESeq OTU differential abundance testing), the OTU matrix from otu picking must be normalized to account for uneven column (sample) sums that are a result of most modern sequencing techniques."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Matrix Normalization""","""For this script, the user supplies an input raw (NOT normalized) OTU matrix (usually always having different column sums), a normalization method (either CSS or DESeqVS), and an output file path, as follows:""","""%prog -i raw_OTU_table.biom -a 'CSS' -o CSS_normalized_OTU_table.biom"""))
script_info['output_description']="""The resulting output file is a normalized count matrix with the same number of OTUs (rows) and samples (columns) as the input raw matrix.  Can be used in all downstream analyses except differential abundance testing, and OTU correlations."""
script_info['required_options']=[\
]
script_info['optional_options']=[
make_option('-i', '--input_path', type='existing_path', help='path to the input raw '
    'matrix file(s) (i.e., the output from OTU picking). Is a directory for '
    'batch processing, filename for a single file operation'),\
make_option('-o', '--out_path', type='new_path', help='output path. directory for '
    'batch processing, filename for single file operation'),\
make_option('-s', '--output_CSS_statistics', default=False, action='store_true', help='output CSS statistics path. '
    'directory for batch processing, filename for single file operation'),\
make_option('-z', '--DESeq_negatives_to_zero', default=False, action='store_true', help='change '
     'negative numbers produced by the DESeq normalization technique'),\
 make_option('-a', '--algorithm', default=False, type='multiple_choice', 
     mchoices=algorithm_list(), help='algorithm to normalize the input raw OTU matrix. '
     'To see the full list of methods and their description use -t'),\
 make_option('-t', '--algorithm_types', action='store_true', default=False, help='show '
     'available OTU matrix normalization algorithms'),
 ]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_path = opts.input_path
    out_path = opts.out_path
    output_CSS_statistics = opts.output_CSS_statistics
    DESeq_negatives_to_zero = opts.DESeq_negatives_to_zero
    algorithm = opts.algorithm
    algorithm_types = opts.algorithm_types


    if algorithm_types:
        print 'Implemented algorithms are:\n%s' % ', '.join(algorithm_list())
    else:
        if len(algorithm)!=1:
            raise ValueError, 'Current implementation only accepts 1 algorithm at a '\
               'time: %s' % algorithm
        else:
            algorithm = algorithm[0]

    if algorithm == 'CSS':            
        if os.path.isdir(input_path):
            multiple_file_normalize_CSS(input_path, out_path, output_CSS_statistics)
        elif os.path.isfile(input_path):
            normalize_CSS(input_path, out_path, output_CSS_statistics)
        else:
            print("io error, check input file path")
        exit(1)


    if algorithm == 'DESeq':            
        if os.path.isdir(input_path):
            multiple_file_normalize_DESeq(input_path, out_path, DESeq_negatives_to_zero)
        elif os.path.isfile(input_path):
            normalize_DESeq(input_path, out_path, DESeq_negatives_to_zero)   
        else:
            print("io error, check input file path")
        exit(1)


if __name__ == "__main__":
    main()



