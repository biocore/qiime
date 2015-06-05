#!/usr/bin/env python
import sys
import os.path
from optparse import OptionParser

from scipy.optimize import curve_fit
from numpy import (exp, cos, pi, tri, argsort, asarray, arange, mean, isnan,
                   zeros, square)

from qiime.parse import parse_distmat

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"


class FitModel(object):

    """This class defines the available models and their functions for a
    semivariogram.
    """

    def __init__(self, x, y, model):
        self.x = x
        self.y = y
        self.model_text = model
        self.model = self._get_model(model)

    # Funcion definition -- defining these in your function makes this
    # very difficult to test
    options = ['nugget', 'exponential', 'gaussian', 'periodic', 'linear']

    def _linear(self, x, a, b):
        self.text = "%f+(%f*x)" % (a, b)
        return a + (b * x)

    def _periodic(self, x, a, b, c):
        self.text = "%f+((1-cos(2*pi*x/%f))*%f)" % (a, b, c)
        return a + ((1 - cos(2 * pi * x / b)) * c)

    def _gaussian(self, x, a, b, c):
        self.text = "%f+((1-exp((-3*x*x)/square(%f)))*%f)" % (a, b, c)
        return a + ((1 - exp((-3 * x * x) / square(b))) * c)

    def _exponential(self, x, a, b, c):
        self.text = "%f+((1-exp(-3*x/%f))*%f)" % (a, b, c)
        return a + ((1 - exp(-3 * x / b)) * c)

    def _nugget(self, x, a):
        self.text = "%f" % (a)
        return a

    def _get_model(self, model):
        if model == 'linear':
            return self._linear
        elif model == 'periodic':
            return self._periodic
        elif model == 'gaussian':
            return self._gaussian
        elif model == 'exponential':
            return self._exponential
        elif model == 'nugget':
            return self._nugget
        else:
            raise ValueError("Unknown model type: %s" % model)

    def __call__(self):
        if self.model_text != 'nugget':
            # what are 3 and 10? should these be parametrizable?
            params, _ = curve_fit(self.model, self.x, self.y)
            y = self.model(self.x, *params)
        else:
            # what are 1 and 1? should these be parametrizable?
            params, _ = curve_fit(self.model, self.x, self.y)
            y = [self.model(self.x, *params)] * len(self.x)

        return y, params, self.text


def hist_bins(bins, vals):
    """ Creates a histogram given the bins and the vals
    :Parameters:
       bins : list
           bins to use
       vals : list
           values to bin

   :Returns:
       bins: array
           The bins
       hist:
           The hist of the values/bins
    """

    hist = zeros(len(bins))
    j = 0
    for i in vals:
        while bins[j] < i:
            j += 1
        hist[j] += 1

    return asarray(bins), hist


def fit_semivariogram(xxx_todo_changeme, xxx_todo_changeme1, model, ranges):
    """ Creates semivariogram values from two distance matrices.
    :Parameters:
       x_file : array matrix distance matrix for x
           distance matrix
       y_file : file handle
           distance matrix file handle
       model: string
           model to fit
       ranges: list
           the list of ranges to bin the data

   :Returns:
       x_vals: array
           Values for x
       y_vals: array
           Values for y
       y_fit: array
           Values for y fitted from model
    """
    (x_samples, x_distmtx) = xxx_todo_changeme
    (y_samples, y_distmtx) = xxx_todo_changeme1
    if x_samples != y_samples:
        lbl_x = list(argsort(x_samples))
        if lbl_x != range(len(lbl_x)):
            tmp = x_distmtx[:, lbl_x]
            x_distmtx = tmp[lbl_x, :]

        lbl_y = list(argsort(y_samples))
        if lbl_y != range(len(lbl_y)):
            tmp = y_distmtx[:, lbl_y]
            y_distmtx = tmp[lbl_y, :]

    # get upper triangle from matrix in a 1d array
    x_tmp_vals = x_distmtx.compress(tri(len(x_distmtx)).ravel() == 0)
    y_tmp_vals = y_distmtx.compress(tri(len(y_distmtx)).ravel() == 0)

    # sorting lists and transforming to arrays
    x_vals, y_vals = [], []
    for i in argsort(x_tmp_vals):
        x_vals.append(x_tmp_vals[i])
        y_vals.append(y_tmp_vals[i])
    x_vals = asarray(x_vals)
    y_vals = asarray(y_vals)

    # fitting model
    fit_func = FitModel(x_vals, y_vals, model)
    y_fit, params, func_text = fit_func()
    x_fit = x_vals

    # section for bins
    if ranges != []:
        # creating bins in x
        min = 0
        x_bins = []
        for r in ranges[:-1]:
            x_bins.extend(arange(min, r[1], r[0]))
            min = r[1]
        x_bins.extend(arange(min, max(x_vals), ranges[-1][0]))
        x_bins[-1] = max(x_vals)

        x_vals, hist = hist_bins(x_bins, x_vals)

        # avg per bin, y values
        y_tmp = []
        for i, val in enumerate(hist):
            if i == 0:
                low = val
                continue
            high = low + val

            y_tmp.append(mean(y_vals[low:high]))

            low = high
        y_vals = asarray(y_tmp)

        # removing nans
        x_vals = x_vals[~isnan(y_vals)]
        y_vals = y_vals[~isnan(y_vals)]

    return x_vals, y_vals, x_fit, y_fit, func_text
