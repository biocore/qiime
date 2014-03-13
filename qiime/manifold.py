#!/usr/bin/env python
import sys
import time
from sklearn import manifold
from numpy import asarray
from biom.parse import parse_biom_table
from biom.table import DenseTable
from qiime.format import format_coords

__author__ = "Joshua Haas"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Joshua Haas,","Antonio Gonzalez Pena","Gregory Ditzler"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Joshua Haas"
__email__ = "laptopdude2@gmail.com"

def compute_manifold(file_name,algorithm):

    otu_table = parse_biom_table(open(file_name,"U"))

    samples = otu_table.sampleIds

    if isinstance(otu_table, DenseTable):
        otumtx = otu_table._data.T
    else:
        otumtx = asarray([v for v in otu_table.iterSampleData()])

    if algorithm=="isomap":
        print otumtx
        fit = manifold.Isomap(n_components=3).fit_transform(otumtx)
    else:
        print("arg in error, unknown algorithm")
        exit(1)

    eigvals = [1,2,3]
    pcnts = [30,20,10]
    
    return format_coords(samples, fit, eigvals, pcnts)

def multiple_file_manifold(input_dir, output_dir, algorithm):
    """perform manifoldss on all distance matrices in the input_dir
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')
                                                        or os.path.isdir(fname))]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        infile = os.path.join(input_dir, fname)
        manifold_res_string = compute_manifold(infile,algorithm)
        outfile = os.path.join(output_dir, algorithm + '_' + base_fname + '.txt')
        outfile = open(outfile, 'w')
        outfile.write(manifold_res_string)
        outfile.close()
