#!/usr/bin/env python
import numpy
import sys
import os.path
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from qiime.format import format_coords
from qiime.parse import parse_distmat

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena", \
    "Catherine Lozupone", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"


def pcoa(file):
    samples, distmtx = parse_distmat(file)
    # coords, each row is an axis
    coords, eigvals = ms.principal_coordinates_analysis(distmtx)
    
    pcnts = (numpy.abs(eigvals)/sum(numpy.abs(eigvals)))*100
    idxs_descending = pcnts.argsort()[::-1]
    coords = coords[idxs_descending]
    eigvals = eigvals[idxs_descending]
    pcnts = pcnts[idxs_descending]
    
    return format_coords(samples, coords.T, eigvals, pcnts)

def multiple_file_pcoa(input_dir, output_dir):
    """perform PCoAs on all distance matrices in the input_dir 
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')\
        or os.path.isdir(fname))]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        infile = os.path.join(input_dir, fname)
        lines = open(infile, 'U')
        pcoa_res_string = pcoa(lines)
        outfile = os.path.join(output_dir, 'pcoa_'+base_fname+'.txt')
        outfile = open(outfile, 'w')
        outfile.write(pcoa_res_string)
        outfile.close()

