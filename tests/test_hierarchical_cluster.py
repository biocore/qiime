#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["justin kuczynski"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from tempfile import mkstemp
from os import close

from qiime.hierarchical_cluster import single_file_upgma, single_file_nj
from qiime.format import format_distance_matrix
from unittest import TestCase, main
import os
import numpy


def remove_files(list_of_filepaths, error_on_missing=True):
    missing = []
    for fp in list_of_filepaths:
        try:
            os.remove(fp)
        except OSError:
            missing.append(fp)

    if error_on_missing and missing:
        raise OSError("Some filepaths were not accessible: %s"
                      % '\t'.join(missing))


class FunctionTests(TestCase):

    """Tests"""

    def setUp(self):
        self._paths_to_clean_up = []

    def test_single_file_upgma(self):
        """ single_file_upgma should throw no errors"""

        titles = ['hi', 'ho']
        distdata = numpy.array([[0, .5], [.5, 0.]])
        fd, fname = mkstemp(prefix='upgma_', suffix='.txt')
        close(fd)
        f = open(fname, 'w')
        self._paths_to_clean_up.append(fname)
        f.write(format_distance_matrix(titles, distdata))
        f.close()

        fd, fname2 = mkstemp(prefix='upgma_', suffix='.txt')
        close(fd)
        self._paths_to_clean_up.append(fname2)
        single_file_upgma(fname, fname2)
        assert(os.path.exists(fname2))

    def test_single_file_nj(self):
        """ single_file_nj should throw no errors"""

        titles = ['hi', 'ho', 'yo']
        distdata = numpy.array([[0, .5, .3], [.5, 0., .9], [.3, .9, 0.]])
        fd, fname = mkstemp(prefix='nj_', suffix='.txt')
        close(fd)
        f = open(fname, 'w')
        self._paths_to_clean_up.append(fname)
        f.write(format_distance_matrix(titles, distdata))
        f.close()

        fd, fname2 = mkstemp(prefix='nj_', suffix='.txt')
        close(fd)
        self._paths_to_clean_up.append(fname2)
        single_file_nj(fname, fname2)
        assert(os.path.exists(fname2))

    def tearDown(self):
        remove_files(self._paths_to_clean_up)


# run tests if called from command line
if __name__ == '__main__':
    main()
