#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["justin kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from qiime.collate_alpha import make_output_row
from unittest import TestCase, main
import os
import numpy


class FunctionTests(TestCase):

    """Tests of top-level functions"""

    def test_make_output_rows(self):
        f_metrics = ['met1']
        metric = 'met1'
        f_samples = ['s1', 's2']
        f_data = numpy.array([[.4], [.8]])
        fname = 'alpha_rarefaction_10_7'
        num_cols = 2
        all_samples = ['s1', 's2']
        res = make_output_row(f_metrics, metric, f_samples,
                              f_data, fname, num_cols, all_samples)
        self.assertEqual(res, ['alpha_rarefaction_10_7', 10, 7, '0.4', '0.8'])

# run tests if called from command line
if __name__ == '__main__':
    main()
