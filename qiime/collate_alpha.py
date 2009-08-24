#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski"] #remember to add yourself
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

from pipe454.parse import (parse_otus, filter_otus_by_lineage, parse_matrix,
    parse_rarefaction_fname)
from pipe454.format import format_otu_table, format_matrix
from pipe454.util import FunctionWithParams

def main(input_dir, output_dir):
    file_names = os.listdir(input_dir)
    example_file = os.path.join(input_dir, file_names[0])
    f = open(example_file, 'r')
    all_metrics, all_samples, example_data = parse_matrix(f)
    num_cols = len(all_samples)
    f.close()

    # make the table 1 row at a time
    # we're building a rarefaction by sample mtx from
    # a sample by metric matrix
    for metric in all_metrics:
        metric_file_data = []
        for fname in file_names:
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
                    print "warning, didn't find sample in example file"

            output_row.insert(1, seqs)
            output_row.insert(2, iter)
            metric_file_data.append(output_row)
        # now have matrix where output_row is rarefaction analysis
        metric_file_data = sorted(metric_file_data, key=operator.itemgetter(1,2))
        row_names = [row.pop(0) for row in metric_file_data]
        col_names = ['sequences per sample', 'iteration'] + all_samples
        out_str = format_matrix(numpy.array(metric_file_data), row_names, \
            col_names)
        f = open(os.path.join(output_dir,metric+'.txt'),'w')
        f.write(out_str)
        f.close()

def make_cmd_parser():
    """returns command-line options"""
    usage = """python collate_alpha.py -i TEST/rare_chao1_PD -o TEST/rare_collated
this creates the files TEST/rare_collated/chao1.txt (and PD_whole_tree.txt)
each is a matrix of rarefaction by sample."""

    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--in_path', dest='input_path', default=None,
        help='output path (a directory)')
    parser.add_option('-o', '--out_path', dest='output_path', default=None,
        help='output path (a directory)')
    options, args = parser.parse_args()
    return options, args


if __name__ == '__main__':
    options,args = make_cmd_parser()
    main(options.input_path, options.output_path)
