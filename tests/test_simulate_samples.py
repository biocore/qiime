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
        clip_on = True
        ints = True
        sparsity = None
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, ints, clip_on, sparsity)
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
            test_std, ints, clip_on, sparsity)
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
        clip_on = False
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, ints, clip_on, sparsity)
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
        clip_on = True
        ints = False
        sparsity = .8
        test_mean = 10.4
        test_std = 1.4
        actual_out = null_from_normal(self.samples, self.otus, test_mean,
            test_std, ints, clip_on, sparsity)
        expected_out = \
            array([[ 10.30446175,  12.79867981,   9.35734325,   9.24298605,
                     10.26216647],
                   [  9.4711304 ,  11.97729029,   8.88809589,   8.79354389,
                      9.78705194],
                   [  9.70275457,  13.10134488,  11.72918913,  10.52257174,
                      8.68439027],
                   [ 11.58210817,   8.99969851,   8.23732046,  12.06324171,
                     10.84371966],
                   [ 11.68920235,  10.84621871,  11.59956286,   9.48856417,
                      8.95206002],
                   [ 11.35423233,   9.27522647,   9.43463031,   9.7622545 ,
                     10.42447082],
                   [  9.90440852,   8.47506819,   9.49893424,   7.28723559,
                     11.27532403],
                   [  8.15711928,   8.85386332,  10.47303111,   9.36461181,
                     12.56022043],
                   [  8.59000033,  10.77387122,  10.34500405,   8.7646691 ,
                     11.13258732],
                   [ 10.15983514,  11.48050677,  11.55290582,  13.42853033,
                     12.27113913]])
        assert_allclose(expected_out, actual_out)
        # using assert_allclose because floats are truncated

    def test_null_from_exponential(self):
        """tests null_from_exponential works properly when seeded"""
        scale = 5
        ints = True
        actual_out = null_from_exponential(self.samples, self.otus, scale, ints)
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
        ints = False
        actual_out = null_from_exponential(self.samples, self.otus, scale, ints)
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