#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

"""runs upgma on a distance matrix file, writes the resulting tree
with distances to specified output file
"""
from optparse import OptionParser
from qiime.parse import parse_distmat
from cogent.core.tree import PhyloNode
from cogent.cluster.UPGMA import UPGMA_cluster
import qiime.hierarchical_cluster
import os.path

def multiple_file_upgma(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    #upgma_script = qiime.hierarchical_cluster.__file__ #removed below
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        single_file_upgma(os.path.join(input_dir, fname),
            os.path.join(output_dir,'upgma_'+base_fname+'.tre'))
            
        # upgma_cmd = 'python ' + upgma_script + ' -i '+\
        #             os.path.join(input_dir, fname) + ' -o ' +\
        #             os.path.join(output_dir,'upgma_'+base_fname+'.tre')
        #         os.system(upgma_cmd)

def single_file_upgma(input_file, output_file):
    # read in dist matrix
    f = open(input_file, 'U')
    headers, data = parse_distmat(f)
    f.close()
    
    # do upgma
    nodes = map(PhyloNode, headers)
    BIG = 1e305
    U = data.copy()
    for i in range(len(U)):
        U[i,i] = BIG
    c = UPGMA_cluster(U, nodes, BIG)

    # write output
    f = open(output_file,'w')
    f.write(c.getNewick(with_distances=True))
    f.close()