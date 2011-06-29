#!/usr/bin/env python


__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Release"

from qiime.plot_semivariogram import hist_bins, fit_semivariogram
from cogent.util.unit_test import TestCase, main
from numpy import asarray

class FunctionTests(TestCase):
    """Tests of top-level functions"""

    def test_hist_bins(self):
        """ test hist_bins """
        x = asarray([3., 4.12310563, 4.24264069, 4.47213595, 5., 5., 5., 5., 5.38516481, 5.65685425, 6.40312424, 6.40312424, 6.70820393, 7.,7.07106781, 7.07106781, 7.28010989, 7.81024968, 8., 8.06225775, 8.06225775, 8.24621125, 9., 9.48683298, 9.48683298, 9.89949494, 9.89949494, 10., 10.04987562, 10.04987562])
        
        bins = [2.0, 5.0, 7.5, 10.0, 11.0]
        hist_res = [0., 8., 9., 11., 2.]
        
        vals, hist = hist_bins(bins , x)
        
        self.assertEqual(vals, bins)
        self.assertEqual(hist, hist_res)
        
        
    def test_fit_semivariogram(self):
        """ test fit_semivariogram """
        x = asarray([[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],[0.0, 0.0, 6.0, 7.0, 8.0, 9.0],[0.0, 0.0, 0.0, 10.0, 11.0, 12.0],[0.0, 0.0, 0.0, 0.0, 13.0, 14.0],[0.0, 0.0, 0.0, 0.0, 0.0, 15.0]])
        y = asarray([[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],[0.0, 0.0, 6.0, 7.0, 8.0, 9.0],[0.0, 0.0, 0.0, 10.0, 11.0, 12.0],[0.0, 0.0, 0.0, 0.0, 13.0, 14.0],[0.0, 0.0, 0.0, 0.0, 0.0, 15.0]])
        ranges = [[2.5, 10.0], [10.0, 400.0], [50.0, -1.0]]
        
        ####
        #### potentially we could create test for each model (TODO?)
        model = "exponential"
        
        # testing without ranges, not checking the fitting cause that is tested somewhere else
        x_vals, y_vals, x_fit, y_fit = fit_semivariogram(x, y, model, [])
        x_vals_exp = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 7.0]
        y_vals_exp = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 7.0]
        self.assertEqual(x_vals, x_vals_exp)
        self.assertEqual(y_vals, y_vals_exp)
        
        # testing with ranges, not checking the fitting cause that is tested somewhere else
        model = "exponential"
        x_vals, y_vals, x_fit, y_fit = fit_semivariogram(x, y, model, ranges)
        x_vals_exp = [0.0, 2.5, 5.0]
        y_vals_exp = [1.5, 3.5, 6.5]
        self.assertEqual(x_vals, x_vals_exp)
        self.assertEqual(y_vals, y_vals_exp)
        

#run tests if called from command line
if __name__ == '__main__':
    main()
    