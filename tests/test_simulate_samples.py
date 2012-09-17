#!/usr/bin/env python
# File created on 11 Sep 2012
from __future__ import division

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"


from cogent.util.unit_test import TestCase, main
from shutil import rmtree
from os.path import exists, join
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, create_dir
from qiime.util import (get_qiime_temp_dir, 
                        get_tmp_filename)
from qiime.test import initiate_timeout, disable_timeout
from qiime.simulate_samples import (null_from_normal, null_from_exponential,
    null_from_data)
from numpy.random import seed
from numpy import array
from numpy.testing import assert_allclose

class TopLevelTests(TestCase):
    
    def setUp(self):
        """setup data used by all functions for testing"""
        self.samples = 5
        self.otus = 10
        seed(0) # this will seed numpy prng at 0 before each test

    def test_null_from_normal(self):
      
        """test that null from normal performs as expected when seeded"""
        test_mean = 10
        test_std = 5
        clip_on = False
        floats = False
        sparsity = None
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, sparsity, floats, clip_on)
        expected_out = \
               array([[ 19.,  12.,  15.,  21.,  19.],
                      [  5.,  15.,   9.,   9.,  12.],
                      [ 11.,  17.,  14.,  11.,  12.],
                      [ 12.,  17.,   9.,  12.,   6.],
                      [  0.,  13.,  14.,   6.,  21.],
                      [  3.,  10.,   9.,  18.,  17.],
                      [ 11.,  12.,   6.,   0.,   8.],
                      [ 11.,  16.,  16.,   8.,   8.],
                      [  5.,   3.,   1.,  20.,   7.],
                      [  8.,   4.,  14.,   2.,   9.]])
        self.assertEqual(expected_out, actual_out)
        # test where values are going to get clipped
        test_mean = 5
        test_std = 5
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, sparsity, floats, clip_on)
        expected_out = \
            array([[  1.,   7.,   2.,   0.,   5.],
                   [  7.,   5.,   7.,   2.,   3.],
                   [  2.,   3.,   1.,   0.,   6.],
                   [  3.,   0.,   7.,   0.,   5.],
                   [  9.,   6.,  11.,   0.,   7.],
                   [  2.,   1.,   2.,   3.,   5.],
                   [  0.,  10.,   7.,   0.,  12.],
                   [ 14.,  11.,   4.,   0.,  10.],
                   [  3.,  11.,   6.,  10.,   7.],
                   [  9.,   5.,  14.,   6.,   7.]])
        self.assertEqual(expected_out, actual_out)
        # test without clipping
        clip_on = True
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, sparsity, floats, clip_on)
        expected_out = \
            array([[ 14.,  -2.,  -1.,  10.,  -1.],
                   [ 15.,   3.,   1.,  15.,  12.],
                   [ 14.,  10.,   1.,  15.,   4.],
                   [  9.,  10.,   4.,   8.,  10.],
                   [  7.,   0.,   6.,  12.,   2.],
                   [  4.,   3.,  14.,   8.,   7.],
                   [  1.,   8.,   2.,   5.,   2.],
                   [  8.,   8.,   4.,   7.,   0.],
                   [ -2.,   7.,   6.,   8.,  17.],
                   [ 10.,   0.,  11.,  -2.,   3.]])
        self.assertEqual(expected_out, actual_out)
        # test with sparsity and ints off
        seed(0)
        clip_on = True
        floats = True
        sparsity = .8
        test_mean = 10.4
        test_std = 1.4
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, sparsity, floats, clip_on)
        expected_out = \
            array([[ 12.86967328,   0.        ,   0.        ,   0.        ,   0.        ],
                   [  9.03181097,   0.        ,   0.        ,   0.        ,   0.        ],
                   [  0.        ,   0.        ,   0.        ,  10.57034502,   0.        ],
                   [ 10.86714406,   0.        ,   0.        ,   0.        ,   0.        ],
                   [  0.        ,  11.31506603,   0.        ,   9.36096897,   0.        ],
                   [  0.        ,   0.        ,   0.        ,   0.        ,   0.        ],
                   [  0.        ,  10.92942753,   0.        ,   0.        ,   0.        ],
                   [ 10.61888856,   0.        ,   0.        ,   0.        ,   0.        ],
                   [  0.        ,   8.41197489,   0.        ,   0.        ,   0.        ],
                   [  9.78669598,   0.        ,   0.        ,   0.        ,   0.        ]])
        assert_allclose(expected_out, actual_out)
        # using assert_allclose because floats are truncated

    def test_null_from_exponential(self):
        """tests null_from_exponential works properly when seeded"""
        scale = 5
        floats = False
        actual_out = null_from_exponential(self.samples, self.otus, scale, floats)
        expected_out = \
            array([[  4.,   6.,   5.,   4.,   3.],
                   [  5.,   3.,  11.,  17.,   2.],
                   [  8.,   4.,   4.,  13.,   0.],
                   [  0.,   0.,   9.,   8.,  10.],
                   [ 19.,   8.,   3.,   8.,   1.],
                   [  5.,   1.,  14.,   4.,   3.],
                   [  2.,   7.,   3.,   4.,   0.],
                   [  5.,   5.,   5.,  14.,   6.],
                   [  2.,   3.,   6.,   0.,   5.],
                   [  6.,   1.,   1.,   2.,   2.]])
        self.assertEqual(expected_out, actual_out)
        #test with ints off
        floats = True
        actual_out = null_from_exponential(self.samples, self.otus, scale, floats)
        expected_out = \
            array([[  4.2221389 ,   2.88662155,  22.2724869 ,   0.53817556,
                      1.17150758],
                   [  0.87956777,   5.29371362,   1.46040267,   3.13970791,
                      1.40138507],
                   [  0.86563726,   0.58477706,   5.34036095,   0.74356136,
                      1.09440301],
                   [  2.30006983,   8.60165826,   0.51072443,   9.09909462,
                      0.50517391],
                   [ 18.74515725,   3.16168301,  18.80963586,   4.64239251,
                      6.72122632],
                   [  0.19988151,   1.66205123,   0.64028381,   1.75588043,
                      0.63194321],
                   [  1.91350479,   2.67442193,   0.33148698,   5.89594764,
                      4.18048773],
                   [  1.5420742 ,   3.70379476,   0.49325157,   4.28947821,
                     13.24627963],
                   [  1.91780105,   5.50422961,   0.70665357,   6.29966916,
                      1.70827085],
                   [  1.01175218,   4.41564523,   0.10156227,   8.82870539,
                      0.02353267]])
        assert_allclose(expected_out, actual_out, rtol=1e-06)
        #using allclose because floats get truncated

    def test_null_from_data(self):
        """tests that null_from_data works as expected when R is seeded"""
        # define prior data that R will use
        data = array([[  4.,   6.,   5.,   4.,   3.],
                      [  5.,   3.,  11.,  17.,   2.],
                      [  8.,   4.,   4.,  13.,   0.],
                      [  0.,   0.,   9.,   8.,  10.],
                      [ 19.,   8.,   3.,   8.,   1.],
                      [  5.,   1.,  14.,   4.,   3.],
                      [  2.,   7.,   3.,   4.,   0.],
                      [  5.,   5.,   5.,  14.,   6.],
                      [  2.,   3.,   6.,   0.,   5.],
                      [  6.,   1.,   1.,   2.,   2.]])
        Rseed = 0
        tpk = 10
        actual_out = null_from_data(data, tpk, Rseed=Rseed)
        expected_out = array([[ 18.,   1.,   0.,   3.,   0.],
                              [  2.,  10.,  10.,   3.,  22.],
                              [  1.,   8.,   4.,   1.,  12.],
                              [ 10.,   1.,   5.,   1.,   2.],
                              [  8.,   0.,  14.,  28.,   5.],
                              [  5.,  16.,   3.,  11.,   0.],
                              [  0.,   2.,   8.,   2.,   0.],
                              [  0.,  12.,   2.,   1.,   3.],
                              [  0.,   0.,   6.,   2.,   5.],
                              [  8.,   2.,   0.,   0.,   3.]])
        self.assertEqual(expected_out, actual_out)


if __name__ == "__main__":
    main()