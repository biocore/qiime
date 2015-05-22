#!/usr/bin/env python
# File created on 27 Sep 2011

from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena", "Andrew J. King"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

"""Contains code for performing distance matrix analyses of a column from a mapping file.

command line usage help: python distance_matrix_from_mapping.py -h

This module has the responsibility for creating a distance matrix from sample to sample using
a column from a mapping file:
    - a tab-delimited mapping file
    - the header of the column to relate

The output is a matrix of distances, incl. row/col headers.
    Note that parser expects first field to be blank, i.e. first char of file
    is expected to be a tab.
"""

from qiime.util import FunctionWithParams
from numpy import array, reshape, nan, radians, zeros
from math import atan, tan, sin, cos, pi, sqrt, atan2, acos, asin


def compute_distance_matrix_from_metadata(column_data):
    """ calculates distance matrix on a single column of a mapping file

    inputs:
     column_data (list of values)
    """
    data_row = array(column_data)
    data_col = reshape(data_row, (1, len(data_row)))
    dist_mtx = abs(data_row - data_col.T)

    return dist_mtx


def dist_vincenty(lat1, lon1, lat2, lon2, iterations=20):
    """Returns distance in meters between two lat long points

       Vincenty's formula is accurate to within 0.5mm, or 0.000015" (!),
       on the ellipsoid being used. Calculations based on a spherical model,
       such as the (much simpler) Haversine, are accurate to around 0.3%
       (which is still good enough for most purposes, of course).
       from: http://www.movable-type.co.uk/scripts/latlong-vincenty.html

       Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the */
       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    */
       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

       This code was modified from geopy and movable-type:
       http://code.google.com/p/geopy/source/browse/trunk/geopy/distance.py?r=105
       http://www.movable-type.co.uk/scripts/latlong-vincenty.html
    """
    if lat1 < -90 or lat1 > 90 or lat2 < -90 or lat2 > 90 or lon1 < -180 or lon1 > 180 or lon2 < -180 or lon2 > 180:
        raise ValueError(
            "Latitude values shoulds range from (-90,90) and longitude from (-180,180) but one of the input values is out of bounds. Latitude_1: %f, Logitude_1: %f, Latitude_2: %f, Logitude_2: %f" %
            (lat1, lon1, lat2, lon2))

    major, minor, f = 6378137, 6356752.314245, 1 / 298.257223563

    lat1, lng1, lat2, lng2 = radians(
        lat1), radians(lon1), radians(lat2), radians(lon2)
    delta_lng = lng2 - lng1
    reduced_lat1, reduced_lat2 = atan(
        (1 - f) * tan(lat1)), atan((1 - f) * tan(lat2))

    sin_reduced1, cos_reduced1 = sin(reduced_lat1), cos(reduced_lat1)
    sin_reduced2, cos_reduced2 = sin(reduced_lat2), cos(reduced_lat2)

    lambda_lng = delta_lng
    lambda_prime = 2 * pi
    while abs(lambda_lng - lambda_prime) > 10e-12 and iterations > 0:
        sin_lambda_lng, cos_lambda_lng = sin(lambda_lng), cos(lambda_lng)

        sin_sigma = sqrt(
            (cos_reduced2 * sin_lambda_lng) ** 2 +
            (cos_reduced1 * sin_reduced2 -
             sin_reduced1 * cos_reduced2 * cos_lambda_lng) ** 2
        )
        if sin_sigma == 0:
            return 0  # Coincident points

        cos_sigma = (
            sin_reduced1 * sin_reduced2 +
            cos_reduced1 * cos_reduced2 * cos_lambda_lng
        )
        sigma = atan2(sin_sigma, cos_sigma)

        sin_alpha = (cos_reduced1 * cos_reduced2 * sin_lambda_lng / sin_sigma)
        cos_sq_alpha = 1 - sin_alpha ** 2

        if cos_sq_alpha != 0:
            cos2_sigma_m = cos_sigma - 2 * \
                (sin_reduced1 * sin_reduced2 / cos_sq_alpha)
        else:
            cos2_sigma_m = 0.0  # Equatorial line

        C = f / 16. * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))

        lambda_prime = lambda_lng
        lambda_lng = (
            delta_lng + (1 - C) * f * sin_alpha * (
                sigma + C * sin_sigma * (
                    cos2_sigma_m + C * cos_sigma * (-1 + 2 * cos2_sigma_m ** 2)
                )
            )
        )
        iterations -= 1

    if iterations == 0:
        raise ValueError("Vincenty formula failed to converge!")

    u_sq = cos_sq_alpha * (major ** 2 - minor ** 2) / minor ** 2
    A = 1 + u_sq / 16384. * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    B = u_sq / 1024. * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    delta_sigma = B * sin_sigma * (
        cos2_sigma_m + B / 4. * (cos_sigma * (-1 + 2 * cos2_sigma_m ** 2) -
                                 B / 6. * cos2_sigma_m * (-3 + 4 * sin_sigma ** 2) *
                                 (-3 + 4 * cos2_sigma_m ** 2))
    )
    s = minor * A * (sigma - delta_sigma)

    return round(s, 3)  # round to 1mm precision


def calculate_dist_vincenty(latitudes, longitudes):
    """Returns the distance matrix from calculating dist_Vicenty

       latitudes, longitudes: list of values, have to be the same size
    """
    assert len(latitudes) == len(
        longitudes), "latitudes and longitudes must be lists of exactly the same size"

    size = len(latitudes)
    dtx_mtx = zeros([size, size])

    for i in range(size):
        for j in range(i, size):
            dtx_mtx[i,
                    j] = dtx_mtx[j,
                                 i] = dist_vincenty(
                latitudes[i],
                longitudes[i],
                latitudes[j],
                longitudes[j])

    return dtx_mtx
