#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division
import operator
import numpy
import os
import sys
from qiime.collate_alpha import write_output_file, make_output_row
from qiime.parse import parse_matrix, parse_rarefaction_fname
from qiime.util import FunctionWithParams
from qiime.util import parse_command_line_parameters, make_option


__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

# collate_alpha.py
script_info = {}
script_info['brief_description'] = """Collate alpha diversity results"""
script_info['script_description'] = """When performing batch analyses on the OTU table (e.g. rarefaction followed by alpha diversity), the result of alpha_diversity.py comprises many files, which need to be concatenated into a single file for generating rarefaction curves.  This script joins those files.
Input files are:
each file represents one (rarefied) otu table
each row in a file represents one sample
each column in a file represents one diversity metric

Output files are:
each file represents one diversity metric
each row in a file represents one (rarefied) otu table
each column in a file represents one sample

The input directory should contain only otu tables. The output directory should be empty or nonexistant and the example file is optional.

If you have a set of rarefied OTU tables, make sure the example file contains every sample present in the otu talbes. You should typically choose the file with the fewest sequences per sample, to avoid files with sparse samples omitted.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """The user inputs the results from batch alpha diversity (e.g. alpha_div/) and the location where the results should be written (e.g. collated_alpha/), as shown by the following command:""",
     """%prog -i alpha_div/ -o collated_alpha/"""))
script_info['output_description'] = """This script takes the resulting files from batch alpha diversity and collates them into (one file for each metric used).

This script transforms a series of files, named (e.g. alpha_rarefaction_20_0.txt, alpha_rarefaction_20_1.txt, etc.) into a (usually much smaller) set of files named (e.g. chao1.txt, PD_whole_tree.txt, etc.), where the columns correspond to samples and rows to the rarefaction files inputted, as shown by the following:

==========================  ====================    =========   ======  ======
\                           sequences per sample    iteration   PC.354  PC.355
==========================  ====================    =========   ======  ======
alpha_rarefaction_20_0.txt  20                      0           0.925   0.915
alpha_rarefaction_20_1.txt  20                      1           0.9     0.89
alpha_rarefaction_20_2.txt  20                      2           0.88    0.915
alpha_rarefaction_20_3.txt  20                      3           0.91    0.93
...                         ...                     ...         ...     ...
==========================  ====================    =========   ======  ======

"""
script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_path',
                help='input path (a directory)'),
    make_option('-o', '--output_path', type='new_dirpath',
                help='output path (a directory).  will be created if needed')
]

script_info['optional_options'] = [
    make_option('-e', '--example_path', type='existing_filepath',
                help='example alpha_diversity analysis file, containing all samples' +
                ' and all metrics to be included in the collated result' +
                '[Default: chosen automatically (see usage string)]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if len(args) != 0:
        parser.error("Positional argument detected.  make sure all" +
                     ' parameters are identified.' +
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
        sorted_fname_table = sorted(
            file_name_table,
            key=operator.itemgetter(1))
        # now map back to file name
        example_fname = file_names[
            file_name_table.index(sorted_fname_table[0])]
        example_filepath = os.path.join(input_dir, example_fname)
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
            f = open(os.path.join(input_dir, fname), 'U')
            f_metrics, f_samples, f_data = parse_matrix(f)
            f.close()
            metric_file_data.append(
                make_output_row(f_metrics, metric, f_samples,
                                f_data, fname, num_cols, all_samples))

        write_output_file(metric_file_data, output_dir, metric, all_samples)


if __name__ == "__main__":
    main()
