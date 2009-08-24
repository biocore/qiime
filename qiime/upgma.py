#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""runs upgma on a distance matrix file, writes the resulting tree
with distances to specified output file
"""
from optparse import OptionParser


from pipe454.parse import parse_distmat

from cogent.core.tree import PhyloNode
from cogent.cluster.UPGMA import UPGMA_cluster
import pipe454.upgma
import os.path

def multiple_file_upgma(options, args):
    upgma_script = pipe454.upgma.__file__
    file_names = os.listdir(options.input_path)
    for fname in file_names:
        upgma_cmd = 'python ' + upgma_script + ' -i '+\
            os.path.join(options.input_path, fname) + ' -o ' +\
            os.path.join(options.output_path,'upgma_'+fname)
        os.system(upgma_cmd)

def single_file_upgma(options, args):
    # read in dist matrix
    f = open(options.input_path, 'r')
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
    f = open(options.output_path,'w')
    f.write(c.getNewick(with_distances=True))
    f.close()


def make_cmd_parser():
    """returns command-line options"""

    usage = """ %prog [options]
use %prog -h for help.

example: 
python upgma.py -i TEST/beta_unweighted_unifrac.txt -o TEST/sample_cluster.tre
creates file TEST/sample_cluster.tre, newick format of upgma clustering based on unifrac of otu_table.

or batch example: 
python upgma.py -i TEST/rare_unifrac -o TEST/rare_unifrac_upgma
processes every file in rare_unifrac, and creates a file "upgma_" + fname
in rare_unifrac_upgma folder
-o is mandatory here
"""
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--input_path', dest='input_path',
        help='input path ')
    parser.add_option('-o', '--output_path', dest='output_path',
        help='output path')
    options, args = parser.parse_args()
    return options, args


if __name__ == '__main__':
    options, args = make_cmd_parser()
    if os.path.isdir(options.input_path):
        multiple_file_upgma(options, args)
    elif os.path.isfile(options.input_path):
        single_file_upgma(options, args)
    else:
        print("io error")
        exit(1)
