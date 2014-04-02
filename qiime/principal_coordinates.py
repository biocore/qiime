#!/usr/bin/env python
import numpy
import os.path

from skbio.core.distance import SymmetricDistanceMatrix
from skbio.maths.stats.ordination import PCoA

from qiime.format import format_coords

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena",
               "Catherine Lozupone", "Emily TerAvest",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


def pcoa(lines):
    dist_mtx = SymmetricDistanceMatrix.from_file(lines)
    pcoa_obj = PCoA(dist_mtx)
    pcoa_results = pcoa_obj.scores()
    # coords, each row is an axis
    coords = pcoa_results.species
    eigvals = pcoa_results.eigvals

    # coords, eigvals = ms.principal_coordinates_analysis(distmtx)

    pcnts = (numpy.abs(eigvals) / sum(numpy.abs(eigvals))) * 100
    idxs_descending = pcnts.argsort()[::-1]
    coords = coords[idxs_descending]
    eigvals = eigvals[idxs_descending]
    pcnts = pcnts[idxs_descending]

    return format_coords(dist_mtx.ids, coords.T, eigvals, pcnts)


def multiple_file_pcoa(input_dir, output_dir):
    """perform PCoAs on all distance matrices in the input_dir
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')
                                                        or os.path.isdir(fname))]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        infile = os.path.join(input_dir, fname)
        lines = open(infile, 'U')
        pcoa_res_string = pcoa(lines)
        outfile = os.path.join(output_dir, 'pcoa_' + base_fname + '.txt')
        outfile = open(outfile, 'w')
        outfile.write(pcoa_res_string)
        outfile.close()
