#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"

"""runs upgma on a distance matrix file, writes the resulting tree
with distances to specified output file
"""
from optparse import OptionParser
from qiime.parse import parse_distmat
from cogent.core.tree import PhyloNode
from cogent.cluster.UPGMA import UPGMA_cluster
import qiime.hierarchical_cluster
import os.path
import sys

def multiple_file_upgma(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    upgma_script = qiime.hierarchical_cluster.__file__
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        upgma_cmd = 'python ' + upgma_script + ' -i '+\
            os.path.join(input_dir, fname) + ' -o ' +\
            os.path.join(output_dir,'upgma_'+base_fname+'.tre')
        os.system(upgma_cmd)

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

usage_str = """usage: %prog [options] {-i INPUT_PATH -o OUTPUT_PATH}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
python %prog -i TEST/beta_unweighted_unifrac.txt -o TEST/sample_cluster.tre
creates file TEST/sample_cluster.tre, newick format of upgma clustering based on
distance matrix in beta_unweighted_unifrac.txt

or batch example: 
python %prog -i TEST/rare_unifrac -o TEST/rare_unifrac_upgma
processes every file in rare_unifrac, and creates a file "upgma_" + 
"base_fname.tre"
in rare_unifrac_upgma folder
-o is mandatory here, created if doesn't exist

description:
relate samples with UPGMA (resulting in a tree), using a distance matrix.
"""
def parse_command_line_parameters():
    """returns command-line options"""

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
    
    parser.add_option('-i', '--input_path',
        help='input path.  directory for batch processing, '+\
         'filename for single file operation [REQUIRED]')
        
    parser.add_option('-o', '--output_path',
        help='output path. directory for batch processing, '+\
         'filename for single file operation [REQUIRED]')

    opts, args = parser.parse_args()
        
    if len(args) != 0:
        parser.error("positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')
         
    required_options = ['input_path','output_path']
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
    return opts, args


if __name__ == '__main__':
    options, args = parse_command_line_parameters()
    if os.path.isdir(options.input_path):
        multiple_file_upgma(options.input_path,options.output_path)
    elif os.path.isfile(options.input_path):
        single_file_upgma(options.input_path, options.output_path)
    else:
        print("io error, check input file path")
        exit(1)
