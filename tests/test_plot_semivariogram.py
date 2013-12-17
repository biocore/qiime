#!/usr/bin/env python


__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

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
        
        self.assertFloatEqual(vals, bins)
        self.assertFloatEqual(hist, hist_res)
        
    def test_reorder_samples(self):
        """ test that regural and irregular order give the same results """
        model = "linear"
        # Test normal order
        x_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        x = asarray([[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],[0.0, 0.0, 6.0, 7.0, 8.0, 9.0],[0.0, 0.0, 0.0, 10.0, 11.0, 12.0],[0.0, 0.0, 0.0, 0.0, 13.0, 14.0],[0.0, 0.0, 0.0, 0.0, 0.0, 15.0]])
        y_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        y = asarray([[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],[0.0, 0.0, 6.0, 7.0, 8.0, 9.0],[0.0, 0.0, 0.0, 10.0, 11.0, 12.0],[0.0, 0.0, 0.0, 0.0, 13.0, 14.0],[0.0, 0.0, 0.0, 0.0, 0.0, 15.0]])
        vals_exp = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 7.0]
        
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (x_lbl,x), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, vals_exp)
        
        # Test altered
        model = "linear"
        # order = [5, 1, 3, 4, 0, 2]
        x_lbl = ['s6', 's2', 's4', 's5', 's1', 's3']
        x = asarray([[0.0,0.0,0.0,0.0,0.0,0.0],[9.0,0.0,7.0,8.0,0.0,6.0],[14.0,0.0,0.0,13.0,0.0,0.0],[15.0,0.0,0.0,0.0,0.0,0.0],[5.0,1.0,3.0,4.0,0.0,2.0],[12.0,0.0,10.0,11.0,0.0,0.0]])
        y_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        y=asarray([[0.0,1.0,2.0,3.0,4.0,5.0],[0.0,0.0,6.0,7.0,8.0,9.0],[0.0,0.0,0.0,10.0,11.0,12.0],[0.0,0.0,0.0,0.0,13.0,14.0],[0.0,0.0,0.0,0.0,0.0,15.0],[0.0,0.0,0.0,0.0,0.0,0.0]])
        vals_exp = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.]
        
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (y_lbl,y), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, vals_exp)
        
    def test_models_semivariograms(self):
        """ test the semivariogram fitting models """
        # All models should return the same x_vals, y_vals, x_fit
        # because we are using the same x
        x_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        x = asarray([[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],[0.0, 0.0, 6.0, 7.0, 8.0, 9.0],[0.0, 0.0, 0.0, 10.0, 11.0, 12.0],[0.0, 0.0, 0.0, 0.0, 13.0, 14.0],[0.0, 0.0, 0.0, 0.0, 0.0, 15.0]])
        vals_exp = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 7.0]
        
        model = "nugget"
        y_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        y = asarray([[0.0, 5.0, 5.0, 5.0, 5.0, 5.0],[0.0, 0.0, 5.0, 5.0, 5.0, 5.0],[0.0, 0.0, 0.0, 5.0, 5.0, 5.0],[0.0, 0.0, 0.0, 0.0, 5.0, 5.0],[0.0, 0.0, 0.0, 0.0, 0.0, 5.0]])
        y_vals_exp = [2.3000000143667378]*(len(x)*2)
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (x_lbl,x), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, y_vals_exp)
        
        model = "exponential"
        y_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        y = asarray([[0.0, 1.0, 22.0, 33.0, 44.0, 55.0],[0.0, 0.0, 66.0, 77.0, 88.0, 99.0],[0.0, 0.0, 0.0, 1010.0, 1111.0, 1212.0],[0.0, 0.0, 0.0, 0.0, 1313.0, 1414.0],[0.0, 0.0, 0.0, 0.0, 0.0, 1515.0]])
        y_vals_exp = [-9.44475054741e-09, -9.44475054741e-09, -9.44475054741e-09, -9.44475054741e-09, 0.999999998426, 2.0000000063, 3.00000001417, 4.00000002204, 5.99999999885, 6.99999998726]
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (x_lbl,x), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, y_vals_exp)
        
        model = "gaussian"
        y_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        y = asarray([[0.0, 1.0, 22.0, 33.0, 44.0, 55.0],[0.0, 0.0, 66.0, 77.0, 88.0, 99.0],[0.0, 0.0, 0.0, 1010.0, 1111.0, 1212.0],[0.0, 0.0, 0.0, 0.0, 1313.0, 1414.0],[0.0, 0.0, 0.0, 0.0, 0.0, 1515.0]])
        y_vals_exp = [0.17373665, 0.17373665, 0.17373665, 0.17373665, 0.54915494, 1.55978856, 2.91608962, 4.28808694, 6.24510109, 6.74689019]
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (x_lbl,x), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, y_vals_exp)
        
        model = "periodic"
        y_lbl = ['s1', 's2', 's3', 's4', 's5', 's6']
        y = asarray([[0.0, 1.0, 22.0, 33.0, 44.0, 55.0],[0.0, 0.0, 66.0, 77.0, 88.0, 99.0],[0.0, 0.0, 0.0, 1010.0, 1111.0, 1212.0],[0.0, 0.0, 0.0, 0.0, 1313.0, 1414.0],[0.0, 0.0, 0.0, 0.0, 0.0, 1515.0]])
        y_vals_exp = [0.23248033, 0.23248033, 0.23248033, 0.23248033, 0.5528678, 1.45081215, 2.74913327, 4.19164973, 6.39844476, 6.72728412]
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (x_lbl,x), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, y_vals_exp)
        
        model = "linear"
        y_lbl = x_lbl
        y = x
        x_vals, y_vals, x_fit, y_fit, func_text = fit_semivariogram((x_lbl,x), (x_lbl,x), model, [])
        self.assertFloatEqual(x_vals, vals_exp)
        self.assertFloatEqual(y_vals, vals_exp)
        self.assertFloatEqual(x_fit, vals_exp)
        self.assertFloatEqual(y_fit, vals_exp)

#run tests if called from command line
if __name__ == '__main__':
    main()
    