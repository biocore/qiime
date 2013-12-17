#!/usr/bin/env python
# File created on 02 May 2012
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "wdwvt1@gmail.com"
 

from cogent.util.unit_test import TestCase, main
from numpy.random import seed
from numpy import array, arange
from qiime.relatedness_library import (reduce_mtx, nri, nti, mpd, mntd,
    random_mpd, random_mntd)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_mpd(self):
        """Test if mean phylogenetic distance (mpd) is calculated correctly.
        Notes:
         The formula used is found Webb et al. 2002. Paper available at:
         http://www.annualreviews.org/doi/pdf/10.1146/annurev.ecolsys.33.010802.150448
        """
        distmat = array([[0.0, 0.4, 2.0, 1.3],
                         [0.4, 0.0, 1.6, 0.2],
                         [2.0, 1.6, 0.0, 1.1],
                         [1.3, 0.2, 1.1, 0.0]])
        self.assertFloatEqual(6.6/6., mpd(distmat))

    def test_mntd(self):
        """Test if mean nearest taxon distance (mntd) is calculated correctly.
        Notes:
         The formula used is found Webb et al. 2002. Paper available at:
         http://www.annualreviews.org/doi/pdf/10.1146/annurev.ecolsys.33.010802.150448
        """
        distmat = array([[0.0, 0.4, 2.0, 1.3],
                         [0.4, 0.0, 1.6, 0.2],
                         [2.0, 1.6, 0.0, 1.1],
                         [1.3, 0.2, 1.1, 0.0]])
        self.assertFloatEqual(1.9/4., mntd(distmat))

    def test_reduce_mtx(self):
        """Test matrix reduction works correctly."""

        distmat = array([[ 0,  1,  2,  3,  4,  5,  6],
                         [ 7,  8,  9, 10, 11, 12, 13],
                         [14, 15, 16, 17, 18, 19, 20],
                         [21, 22, 23, 24, 25, 26, 27],
                         [28, 29, 30, 31, 32, 33, 34],
                         [35, 36, 37, 38, 39, 40, 41],
                         [42, 43, 44, 45, 46, 47, 48]])
        inds_1 = [4,5,1]
        exp_1 = array([[32, 33, 29],
                       [39, 40, 36],
                       [11, 12,  8]])
        inds_2 = [1,3,6]
        exp_2 = array([[ 8, 10, 13],
                       [22, 24, 27],
                       [43, 45, 48]])
        self.assertEqual(exp_1, reduce_mtx(distmat, inds_1))
        self.assertEqual(exp_2, reduce_mtx(distmat, inds_2))

    def test_random_mpd(self):
        """Test random mpd draws are calculated correctly when seeded."""
        seed(0)

        distmat = array([[ 0.  ,  0.85,  0.5 ,  0.14,  0.36],
                         [ 0.85,  0.  ,  0.79,  0.25,  0.47],
                         [ 0.5 ,  0.79,  0.  ,  0.24,  0.46],
                         [ 0.14,  0.25,  0.24,  0.  ,  0.8 ],
                         [ 0.36,  0.47,  0.46,  0.8 ,  0.  ]])

        # test calculated by hand to ensure correct. first 3 random index draws
        # are: [2,0,1], [2,1,0], [1,4,3], avgs = 2.14/3, 2.14/3, 1.52/3
        # mean= .62, std means =0.13199326582148888
        obs_mean, obs_std = random_mpd(distmat, n=3, iters=3)
        self.assertFloatEqual(obs_mean, 0.64444444444444449)
        self.assertFloatEqual(obs_std, 0.097423600963479878)

    def test_random_mntd(self):
        """Test random mntd draws are calculated correctly when seeded."""
        seed(0)

        distmat = array([[ 0.  ,  0.85,  0.5 ,  0.14,  0.36],
                         [ 0.85,  0.  ,  0.79,  0.25,  0.47],
                         [ 0.5 ,  0.79,  0.  ,  0.24,  0.46],
                         [ 0.14,  0.25,  0.24,  0.  ,  0.8 ],
                         [ 0.36,  0.47,  0.46,  0.8 ,  0.  ]])

        # test calculated by hand to ensure correct. first 3 random index draws
        # are: [2,0,1], [2,1,0], [1,4,3], avgs = 0.5966666666666667, 
        # 0.5966666666666667, .97/3
        obs_mean, obs_std = random_mntd(distmat, n=3, iters=3)
        self.assertFloatEqual(obs_mean, 0.50555555555555554)
        self.assertFloatEqual(obs_std, 0.12885056901621536)

    def test_nri(self):
        """Test that nri works correctly."""
        # using the input of the distance matrix generated from this tree
        # '(((sp1:.06,sp2:.1)A:.031,(sp3:.001,sp4:.01)B:.2)AB:.4,((sp5:.03,sp6:.02)C:.13,(sp7:.01,sp8:.005)D:.1)CD:.3)root;'
        # using this tree and this clumping we get almost exact agreement 
        # between phylocom and relatedness.py
        distmat = array([[ 0.   ,  0.16 ,  0.292,  0.301,  0.951,  0.941,  0.901,  0.896],
                       [ 0.16 ,  0.   ,  0.332,  0.341,  0.991,  0.981,  0.941,  0.936],
                       [ 0.292,  0.332,  0.   ,  0.011,  1.061,  1.051,  1.011,  1.006],
                       [ 0.301,  0.341,  0.011,  0.   ,  1.07 ,  1.06 ,  1.02 ,  1.015],
                       [ 0.951,  0.991,  1.061,  1.07 ,  0.   ,  0.05 ,  0.27 ,  0.265],
                       [ 0.941,  0.981,  1.051,  1.06 ,  0.05 ,  0.   ,  0.26 ,  0.255],
                       [ 0.901,  0.941,  1.011,  1.02 ,  0.27 ,  0.26 ,  0.   ,  0.015],
                       [ 0.896,  0.936,  1.006,  1.015,  0.265,  0.255,  0.015,  0.   ]])

        seed(0)
        obs_nri = nri(distmat, list(arange(8)), [0,1,3,6],1000)
        self.assertFloatEqual(0.45181900877903675, obs_nri)

    def test_nti(self):
        """Test that nti works correctly."""
        # using the input of the distance matrix generated from this tree
        # '(((sp1:.06,sp2:.1)A:.031,(sp3:.001,sp4:.01)B:.2)AB:.4,((sp5:.03,sp6:.02)C:.13,(sp7:.01,sp8:.005)D:.1)CD:.3)root;'
        # using this tree and this clumping we get almost exact agreement 
        # between phylocom and relatedness.py
        distmat = array([[ 0.   ,  0.16 ,  0.292,  0.301,  0.951,  0.941,  0.901,  0.896],
                       [ 0.16 ,  0.   ,  0.332,  0.341,  0.991,  0.981,  0.941,  0.936],
                       [ 0.292,  0.332,  0.   ,  0.011,  1.061,  1.051,  1.011,  1.006],
                       [ 0.301,  0.341,  0.011,  0.   ,  1.07 ,  1.06 ,  1.02 ,  1.015],
                       [ 0.951,  0.991,  1.061,  1.07 ,  0.   ,  0.05 ,  0.27 ,  0.265],
                       [ 0.941,  0.981,  1.051,  1.06 ,  0.05 ,  0.   ,  0.26 ,  0.255],
                       [ 0.901,  0.941,  1.011,  1.02 ,  0.27 ,  0.26 ,  0.   ,  0.015],
                       [ 0.896,  0.936,  1.006,  1.015,  0.265,  0.255,  0.015,  0.   ]])

        seed(0)
        obs_nti = nti(distmat, list(arange(8)), [0,1,3,6],1000)
        self.assertFloatEqual(-1.2046544711672049, obs_nti)


#run unit tests if run from command-line
if __name__ == '__main__':
    main()

