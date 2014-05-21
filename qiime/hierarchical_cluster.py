#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""runs upgma or nj on a distance matrix file, writes the resulting tree
with distances to specified output file
"""
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
from optparse import OptionParser
from qiime.parse import parse_distmat
from scipy.cluster.hierarchy import linkage
from skbio.core.tree import TreeNode
from skbio.core.distance import DistanceMatrix
from cogent.phylo.nj import nj
import os.path


def multiple_file_upgma(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        single_file_upgma(os.path.join(input_dir, fname),
                          os.path.join(output_dir, 'upgma_' + base_fname + '.tre'))


def single_file_upgma(input_file, output_file):
    # read in dist matrix
    dist_mat = DistanceMatrix.from_file(input_file)

    # SciPy uses average as UPGMA:
    # http://docs.scipy.org/doc/scipy/reference/generated/
    #    scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
    linkage_matrix = linkage(dist_mat.condensed_form(), method='average')

    tree = TreeNode.from_linkage_matrix(linkage_matrix, dist_mat.ids)

    # write output
    f = open(output_file, 'w')
    try:
        f.write(tree.to_newick(with_distances=True))
    except AttributeError:
        if c is None:
            raise RuntimeError("""input file %s did not make a UPGMA tree.
 Ensure it has more than one sample present""" % (str(input_file),))
        raise
    f.close()


def single_file_nj(input_file, output_file):
    # read in dist matrix
    f = open(input_file, 'U')
    headers, data = parse_distmat(f)
    f.close()

    # do nj
    distdict = {}
    for i in range(len(headers)):
        for j in range(len(headers)):
            distdict[(headers[i], headers[j])] = data[i, j]  # need j,i too?

    tree = nj(distdict)

    # write output
    f = open(output_file, 'w')
    f.write(tree.getNewick(with_distances=True))
    f.close()


def multiple_file_nj(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        single_file_nj(os.path.join(input_dir, fname),
                       os.path.join(output_dir, 'nj_' + base_fname + '.tre'))
