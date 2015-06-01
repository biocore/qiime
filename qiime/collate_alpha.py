#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""takes as input alpha diveristy files, format:
file: one rarefaction depth/iteration
rows: samples
cols: alpha diversity metrics

and collates them into:
file: one alpha metric
rows: rarefaction depth/iteration
cols: samples

example:
python collate_alpha.py -i rare_chao1_PD -o rare_collated
doesn't handle written n/a values well, can't read in
"""
import operator
from optparse import OptionParser
import numpy
import os
import sys

from qiime.parse import parse_matrix, parse_rarefaction_fname
from qiime.format import format_matrix
from qiime.util import FunctionWithParams


def write_output_file(metric_file_data, output_dir, metric, all_samples):
    # now have matrix where output_row is rarefaction analysis
    metric_file_data = sorted(metric_file_data, key=operator.itemgetter(1, 2))
    row_names = [row.pop(0) for row in metric_file_data]
    col_names = ['sequences per sample', 'iteration'] + all_samples
    # Numpy shows weird behaviour when converting metric_file_data to array
    # it truncates some values, so better go with straight list of lists
    # format_matrix() now takes 2d lists as well as arrays.
    out_str = format_matrix(metric_file_data, row_names, col_names)
    f = open(os.path.join(output_dir, metric + '.txt'), 'w')
    f.write(out_str)
    f.close()


def make_output_row(f_metrics, metric, f_samples, f_data,
                    fname, num_cols, all_samples):
    f_col = f_metrics.index(metric)

    # first 3 cols are fname, seqs/sample, iteration
    try:
        base, seqs, iter, ext = parse_rarefaction_fname(fname)
    except:
        seqs, iter = 'n/a', 'n/a'
    output_row = [fname] + ['n/a'] * num_cols
    for f_row, sample in enumerate(f_samples):
        try:
            output_row[all_samples.index(sample) + 1] = \
                str(f_data[f_row, f_col])
        except IndexError:
            print("warning, didn't find sample in example file." +
                  "exiting", sample, fname, metric)
            raise  # re-raise error

    output_row.insert(1, seqs)
    output_row.insert(2, iter)
    return output_row
