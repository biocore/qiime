#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
""" computes cluster quality the default way, bet/within"""


def clust_qual_ratio(dists, map_data, category):
    """ measures within category and between category dissimilarities

    input:
    distance matrix [header, numpy 2d array],
    [map_dict, comments]
    category in mapping data (string)

    returns [dissims between clusters], [dissims within same cluster]
    """
    map_dict, comments = map_data
    dist_headers, dmtx = dists
    if len(set(dist_headers)) != len(dist_headers):
        raise RuntimeError("Error: distance matrix headers are non unique")

    cats = sorted([map_dict[sam][category] for sam in map_dict])
    diff_dists = []
    same_dists = []
    for i in range(len(dist_headers)):
        for j in range(i):
            if map_dict[dist_headers[i]][category] == \
                    map_dict[dist_headers[j]][category]:
                same_dists.append(dmtx[i, j])
            else:
                diff_dists.append(dmtx[i, j])

    return diff_dists, same_dists
