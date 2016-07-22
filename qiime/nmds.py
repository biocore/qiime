#!/usr/bin/env python
from __future__ import division
import numpy
import os.path
import cogent.cluster.nmds as nmds_module
from qiime.format import format_nmds_coords
from qiime.parse import parse_distmat

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


def nmds(file, dimensions=2):
    samples, distmtx = parse_distmat(file)
    nmds_res = nmds_module.NMDS(distmtx, verbosity=0, dimension=dimensions)
    pts = nmds_res.getPoints()
    stress = nmds_res.getStress()

    return format_nmds_coords(samples, pts, stress)


def multiple_file_nmds(input_dir, output_dir, dimensions=2):
    """perform PCoAs on all distance matrices in the input_dir
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        infile = os.path.join(input_dir, fname)
        lines = open(infile, 'U')
        nmds_res_string = nmds(lines, dimensions)
        outfile = os.path.join(output_dir, 'nmds_' + base_fname + '.txt')
        outfile = open(outfile, 'w')
        outfile.write(nmds_res_string)
        outfile.close()
