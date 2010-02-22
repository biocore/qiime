#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division
import operator
import numpy
import os
import sys
from optparse import OptionParser,make_option
from qiime.collate_alpha import write_output_file, make_output_row
from qiime.parse import (parse_otus, filter_otus_by_lineage, parse_matrix,
    parse_rarefaction_fname)
from qiime.format import format_otu_table, format_matrix
from qiime.util import FunctionWithParams
from qiime.util import parse_command_line_parameters


__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 
script_description = """collate_alpha collates output files from alpha_diversity.py"""

script_usage = """usage: %prog [options] {-i INPUT_PATH -o OUTPUT_PATH}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
python %prog -i TEST/rare_chao1_PD -o TEST/rare_collated -e
TEST/rare_chao1_PD/alpha_rarefaction_200_0.txt 
this creates the files TEST/rare_collated/chao1.txt (and PD_whole_tree.txt)
each is a matrix of rarefaction by sample.

input directory should have only otu tables
output dir should be empty or nonexistant
example file is optional.  

If you have a set of rarefied OTU tables,
make sure the example file contains every sample present in the otu talbes.
typically choose the file with the fewest sequences per sample, to avoid
files with sparse samples omitted.  This is the default behavior.
"""

required_options = [\
 make_option('-i', '--input_path',
 help='input path (a directory)'),
 make_option('-o', '--output_path',
 help='output path (a directory).  will be created if needed')
]

optional_options = [\
 make_option('-e', '--example_path',
 help='example alpha_diversity analysis file, containing all samples'+\
 ' and all metrics to be included in the collated result'+\
 '[Default: chosen automatically (see usage string)]')
]




def main():
    option_parser, opts, args = parse_command_line_parameters(\
      script_description=script_description,\
      script_usage=script_usage,\
      version=__version__,\
      required_options=required_options,\
      optional_options=optional_options)
    
    if len(args) != 0:
        parser.error("Positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')

    input_dir = opts.input_path
    output_dir = opts.output_path
    example_filepath = opts.example_path 
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    if example_filepath is None:
        # table row is base_name, seqs_per_sam, iters, ext
        file_name_table = map(parse_rarefaction_fname, file_names)
        # sort on seqs/sam
        sorted_fname_table = sorted(file_name_table, key=operator.itemgetter(1))
        # now map back to file name
        example_fname = file_names[file_name_table.index(sorted_fname_table[0])]
        example_filepath = os.path.join(input_dir,example_fname)
    f = open(example_filepath, 'U')
    all_metrics, all_samples, example_data = parse_matrix(f)
    num_cols = len(all_samples)
    f.close()
    
    # make the table 1 row at a time
    # we're building a rarefaction by sample mtx from
    # a sample by metric matrix
    # each metric is one output file
    for metric in all_metrics:
        metric_file_data = []
        for fname in file_names:
            # f_ here refers to the input file currently being processed
            # to distinguish from the output file we're building
            f = open(os.path.join(input_dir,fname), 'U')
            f_metrics, f_samples, f_data = parse_matrix(f)
            f.close()
            metric_file_data.append(\
              make_output_row(f_metrics, metric, f_samples, 
                f_data, fname,num_cols,all_samples))

        write_output_file(metric_file_data, output_dir, metric, all_samples)

    
if __name__ == "__main__":
    main()
