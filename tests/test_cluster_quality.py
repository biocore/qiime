#!/usr/bin/env python

"""Tests of code for cluster quality"""

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

import numpy
from os import remove

from cogent.util.misc import create_dir
from unittest import TestCase, main
from qiime.util import get_tmp_filename

from qiime.cluster_quality import clust_qual_ratio


class FunctionTests(TestCase):

    """Tests of the abstract OtuPicker class"""

    def test_clust_qual_ratio(self):
        dist_labels = ['s1', 's2', 's3']
        dmtx = numpy.array([[0, 2.1, 3], [2.1, 0, 1], [3, 1, 0]])
        dists = [dist_labels, dmtx]
        map_data = \
            {'s1': {'color': 'red'},
             's2': {'color': 'blue'},
                's3': {'color': 'red'}}
        bet, within = clust_qual_ratio(dists, [map_data, []], 'color')
        self.assertFloatEqual(sorted(bet), [1, 2.1])
        self.assertFloatEqual(within, [3.])

# run unit tests if run from command-line
if __name__ == '__main__':
    main()
