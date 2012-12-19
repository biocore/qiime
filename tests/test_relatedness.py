#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren"]
__license__ = "GPL"
__version__ = "1.6.0"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Release"
 

from cogent.util.unit_test import TestCase, main
from numpy.random import seed
from numpy import array
from qiime.relatedness_library import (mntd, mntd_mean_sd, mpd, mpd_mean_sd, nri, nti,
    take_distmat_data, take_random_ids)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    # nti tests
    def test_nti(self):
        "test that nti returns the correct value when shuffle seeded at 0"
        # define nec components
        seed(0)
        datamtx = array([[0.0, 0.4, 2.0, 1.3],
                         [0.4, 0.0, 1.6, 0.2],
                         [2.0, 1.6, 0.0, 1.1],
                         [1.3, 0.2, 1.1, 0.0]])
        all_ids = [0,1,2,3]
        ids_to_keep = [0,1,2]
        iters = 3

        expected_out = 0.08151391459392235
        actual_out = nti(datamtx, all_ids, ids_to_keep, iters)

        self.assertFloatEqual(expected_out, actual_out)
       
        # test that an unsorted list of ids returns same result
        unsorted_ids_to_keep = [2,0,1]
        seed(0)
        actual_out = nti(datamtx, all_ids, ids_to_keep, iters)

        self.assertFloatEqual(expected_out, actual_out)


        # test for situations where standard deviation would be close to zero
        # but not actually zero due to floating point arithmetic errors
        datamtx = array([[0.0, 0.1, 0.1, 0.1],
                         [0.1, 0.0, 0.1, 0.1],
                         [0.1, 0.1, 0.0, 0.1],
                         [0.1, 0.1, 0.1, 0.0]])
        all_ids = [0,1,2,3]
        ids_to_keep = [0,1,2,3]

        self.assertRaises(ValueError, nti, datamtx, all_ids, ids_to_keep, iters)

    def test_mntd(self):
        """test mntd walks through input matrices correctly"""
        # define datamtx
        datamtx = array([[0.0, 0.4, 2.0, 1.3],
                         [0.4, 0.0, 1.6, 0.2],
                         [2.0, 1.6, 0.0, 1.1],
                         [1.3, 0.2, 1.1, 0.0]])

        expected_out_m = 1.9/4.
        actual_out_m = mntd(datamtx)
        self.assertFloatEqual(expected_out_m, actual_out_m)

    def test_mntd_mean_sd(self):
        """tests that randomizations occur correctly given numpy seed 0"""
        # define data
        datamtx = array([[0.0, 0.4, 2.0, 1.3],
                         [0.4, 0.0, 1.6, 0.2],
                         [2.0, 1.6, 0.0, 1.1],
                         [1.3, 0.2, 1.1, 0.0]])
        seed(0)
        all_ids = [0,1,2,3]
        num_to_take = 3
        iters = 3

        means = [1.5/3., 2.4/3., 3.5/3.]
        expected_out_m = 0.8222222222222223
        expected_out_s = 0.27261875880856218

        actual_out_m, actual_out_s = \
            mntd_mean_sd(datamtx, all_ids, num_to_take, iters)

        self.assertFloatEqual(expected_out_m, actual_out_m)
        self.assertFloatEqual(expected_out_s, actual_out_s)

    # nri tests

    def test_nri(self):
        """tests nri returns correct result with sorted/unsorted taxa list,seed0"""
        # define nec components
        seed(0)
        all_ids = [0,1,2,3,4]
        group_ids = [0,2,3,4]
        iters = 3
        datamtx = array([[0., 19., 2., 13., 1.],
                         [19., 0., 15., 8., 9.],
                         [2., 15., 0., 10., 7.],
                         [13., 8., 10., 0., 12.],
                         [1., 9., 7., 12., 0.]])
        expected_out = 2.70454608694765

        actual_out = nri(datamtx, all_ids, group_ids, iters)

        self.assertFloatEqual(expected_out, actual_out)
       
        # test that an unsorted list of ids returns same result
        unsorted_group_ids = [2,0,4,3]
        seed(0)
        actual_out = nri(datamtx, all_ids, unsorted_group_ids, iters)

        self.assertFloatEqual(expected_out, actual_out)

        # test for situations where standard deviation would be close to zero
        # but not actually zero due to floating point arithmetic errors
        datamtx = array([[0.0, 0.1, 0.1, 0.1],
                         [0.1, 0.0, 0.1, 0.1],
                         [0.1, 0.1, 0.0, 0.1],
                         [0.1, 0.1, 0.1, 0.0]])
        all_ids = [0,1,2,3]
        ids_to_keep = [0,1,2,3]

        self.assertRaises(ValueError, nri, datamtx, all_ids, ids_to_keep, iters)



    def test_mpd_mean_sd(self):
        """tests randomizations occur correctly given numpy seed 0"""
        #define some data
        seed(0)
        all_ids = [0,1,2,3,4]
        num_to_take = 4
        iters = 3
        datamtx = array([[0., 19., 2., 13., 1.],
                         [19., 0., 15., 8., 9.],
                         [2., 15., 0., 10., 7.],
                         [13., 8., 10., 0., 12.],
                         [1., 9., 7., 12., 0.]])
        
        expected_out_m = 10.11111111111111
        expected_out_s = 0.9654526220545977
        actual_out_m, actual_out_s = mpd_mean_sd(datamtx,all_ids,num_to_take,iters)

        self.assertFloatEqual(expected_out_m, actual_out_m)
        self.assertFloatEqual(expected_out_s, actual_out_s)

    def test_mpd(self):
        """tests mpd walks through input matrices correctly"""
        #define data
        dist_mat = array([[ 0,  1,  2,  3,  4],
                          [ 5,  6,  7,  8,  9],
                          [10, 11, 12, 13, 14],
                          [15, 16, 17, 18, 19],
                          [20, 21, 22, 23, 24]])

        expected_out_m = 8
        actual_out_m = mpd(dist_mat)

        self.assertEqual(expected_out_m, actual_out_m)

    # joint tests

    def test_take_random_ids(self):
        """tests shuffle works as expected with seed 0"""
        #define data
        from numpy.random import seed
        seed(0)
        all_ids = [0,1,2,3,4,5,6,7]
        num_to_take = 5
        expected_out = [6,2,1,7,3]
        actual_out = take_random_ids(all_ids, num_to_take)

        self.assertEqual(expected_out, actual_out)
        self.assertRaises(ValueError,
            take_random_ids, all_ids, 10)

    def test_take_distmat_data(self):
        """tests array reduction to just row/cols that are wanted works"""
        # define some data
        dist_mat = array([[ 0,  1,  2,  3,  4],
                          [ 5,  6,  7,  8,  9],
                          [10, 11, 12, 13, 14],
                          [15, 16, 17, 18, 19],
                          [20, 21, 22, 23, 24]])
        all_ids = [1,2,3,4,5]
        ids_to_keep = [1,3,5]
        expected_out = array([[0,2,4],
                              [10,12,14],
                              [20,22,24]])
        actual_out = take_distmat_data(dist_mat,all_ids,ids_to_keep)

        self.assertEqual(expected_out, actual_out)

#run unit tests if run from command-line
if __name__ == '__main__':
    main()

