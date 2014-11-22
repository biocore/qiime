#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

from qiime.util import parse_command_line_parameters, make_option
from qiime.differential_abundance import DA_fitZIG, multiple_file_DA_fitZIG, DA_DESeq2, multiple_file_DA_DESeq2, algorithm_list

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
script_info['brief_description'] = "OTU differential abundance between two sample categories"
script_info['script_description'] = "OTU differential abundance testing is commonly used to identify OTUs that differ between two mapping file sample categories (i.e. L_palm and Tongue body sites)."
# Members of the tuple in script_usage are (title, description, example call)
script_info['script_usage'] = [("""OTU Differential Abundance Testing""","""For this script, the user supplies an input raw (NOT normalized) OTU matrix (usually always having different column sums), an output file, a mapping file, a category in the mapping file for which it is desired to test differential abundance, and the algorithm that the user want for differential abundance testing, as follows:""","""differential_abundance.py -i raw_OTU_table.biom -o DE_OTUs.txt -m mapping_fp.txt -a 'metagenomeSeq_fitZIG' -c 'BODY_SITE' -x 'Tongue' -y 'L_palm'""")]
script_info['output_description']= "The resulting output OTU txt file contains a list of all the OTUs in the input matrix, along with their associated statistics and FDR p-values."
script_info['required_options']=[\
 make_option('-i', '--input_path',type='existing_path', help='path to the input raw '
    'matrix file(s) (i.e., the output from OTU picking). '
    'processing and a filename for a single file operation.'),\
 make_option('-o', '--out_path',type='new_path', help='output path. directory for '
    'batch processing, filename for single file operation'),\
 make_option('-a', '--algorithm', type='multiple_choice', 
     mchoices=algorithm_list(), help='algorithm to calculate the differentially abundant OTUs '
     'To see the full list of methods and their description use -s'),\
 make_option('-m', '--mapping_file_path', type='new_path', help='mapping file path.'),\
 make_option('-c', '--mapping_file_category', type='string', help='mapping file category, e.g. "BODY_SITE".'),\
 make_option('-x', '--mapping_file_subcategory_1', type='string', help='mapping file subcategory, e.g. "L_palm".'),\
 make_option('-y', '--mapping_file_subcategory_2', type='string', help='mapping file subcategory, e.g. "Tongue".'),\
]
script_info['optional_options']=[\
 make_option('-s', '--show_algorithms',action='store_true', default=False, help='show '
     'available OTU differential abundance testing algorithms'),\
 make_option('-d', '--DESeq2_diagnostic_plots', default=False, action='store_true', help='show '
     'a MA plot - y axis: log2 fold change, x axis: average size factor normalized OTU value. '
     'also show a Dispersion Estimate plot - visualize the fitted dispersion vs. mean relationship'),
 ]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_path = opts.input_path
    out_path = opts.out_path
    algorithm = opts.algorithm
    mapping_fp = opts.mapping_file_path
    mapping_category = opts.mapping_file_category
    subcategory_1 = opts.mapping_file_subcategory_1
    subcategory_2 = opts.mapping_file_subcategory_2
    show_algorithms = opts.show_algorithms
    DESeq2_diagnostic_plots = opts.DESeq2_diagnostic_plots


    if show_algorithms:
        print 'Implemented algorithms are:\n%s' % ', '.join(algorithm_list())
    else:
        if len(algorithm)!=1:
            raise ValueError, 'Current implementation only accepts 1 algorithm at a '\
               'time: %s' % algorithm
        else:
            algorithm = algorithm[0]

    if algorithm == 'metagenomeSeq_fitZIG':            
        if os.path.isdir(input_path):
            multiple_file_DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2)   
        elif os.path.isfile(input_path):
            DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2)   
        else:
            print("io error, check input file path")
        exit(1)

    if algorithm == 'DESeq2_nbinom':            
        if os.path.isdir(input_path):
            multiple_file_DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots)   
        elif os.path.isfile(input_path):
            DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots)   
        else:
            print("io error, check input file path")
        exit(1)


if __name__ == "__main__":
    main()
