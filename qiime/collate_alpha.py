#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""takes as input alpha diveristy files, format:
file: one rarefaction depth/iteration
rows: samples
cols: alpha diversity metrics

and collates them into:
file: one alpha metric
rows: rarefaction depth/iteration
cols: samples

example:
python collate_alpha.py -i TEST/rare_chao1_PD -o TEST/rare_collated
doesn't handle written n/a values well, can't read in
"""
import operator
from optparse import OptionParser
import numpy
import os
import sys

from qiime.parse import (parse_otus, filter_otus_by_lineage, parse_matrix,
    parse_rarefaction_fname)
from qiime.format import format_otu_table, format_matrix
from qiime.util import FunctionWithParams

def main(input_dir, output_dir, example_filepath=None):
    if not os.path.exists(options.output_path):
        os.makedirs(options.output_path)
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
    f = open(example_filepath, 'r')
    all_metrics, all_samples, example_data = parse_matrix(f)
    num_cols = len(all_samples)
    f.close()

    # make the table 1 row at a time
    # we're building a rarefaction by sample mtx from
    # a sample by metric matrix
    for metric in all_metrics:
        metric_file_data = []
        for fname in file_names:
            # f_ here refers to the input file currently being processed
            # to distinguish from the output file we're building
            f = open(os.path.join(input_dir,fname), 'r')
            f_metrics, f_samples, f_data = \
                parse_matrix(f)
            f.close()
            f_col = f_metrics.index(metric)

            # first 3 cols are fname, seqs/sample, iteration
            try:
                base, seqs, iter, ext = parse_rarefaction_fname(fname)
            except:
                seqs, iter = 'n/a', 'n/a'
            output_row = [fname] + ['n/a']*num_cols
            for f_row, sample in enumerate(f_samples):
                try:
                    output_row[all_samples.index(sample)+1] = \
                        str(f_data[f_row,f_col])
                except IndexError:
                    print("warning, didn't find sample in example file."+\
                     "exiting", sample, fname, metric)
                    raise # re-raise error
                    

            output_row.insert(1, seqs)
            output_row.insert(2, iter)
            metric_file_data.append(output_row)
        # now have matrix where output_row is rarefaction analysis
        metric_file_data = sorted(metric_file_data,key=operator.itemgetter(1,2))
        row_names = [row.pop(0) for row in metric_file_data]
        col_names = ['sequences per sample', 'iteration'] + all_samples
        out_str = format_matrix(numpy.array(metric_file_data), row_names, \
            col_names)
        f = open(os.path.join(output_dir,metric+'.txt'),'w')
        f.write(out_str)
        f.close()

usage_str = """usage: %prog [options] {-i INPUT_PATH -o OUTPUT_PATH}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
python %prog -i TEST/rare_chao1_PD -o TEST/rare_collated -e TEST/rare_chao1_PD/alpha_rarefaction_200_0.txt 
this creates the files TEST/rare_collated/chao1.txt (and PD_whole_tree.txt)
each is a matrix of rarefaction by sample.

input directory should have only otu tables
output dir should be empty or nonexistant
example file is optional.  if you have a set of rarefied otu talbes,
make sure the example file contains every sample present in the otu talbes.
typically choose the file with the fewest sequences per sample, to avoid
files with sparse samples omitted.  This is the default behavior.
"""
def parse_command_line_parameters():
    """returns command-line options"""

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
    
    parser.add_option('-i', '--input_path',
        help='input path (a directory) [REQUIRED]')
        
    parser.add_option('-o', '--output_path',
        help='output path (a directory).  will be created if needed [REQUIRED]')

    parser.add_option('-e', '--example_path', 
        default=None,
        help='example alpha_diversity analysis file, containing all samples'+\
        ' and all metrics to be included in the collated result'+\
        '[default: chosen automatically (see usage string)]')  

    opts, args = parser.parse_args()
        
    if len(args) != 0:
        parser.error("positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')
         
    required_options = ['input_path','output_path']
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
    return opts, args

if __name__ == '__main__':
    options,args = parse_command_line_parameters()
    main(options.input_path, options.output_path, options.example_path)
