#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

import numpy
import sys
import os
from qiime.parse import parse_distmat
from qiime.format import format_distance_matrix
from optparse import OptionParser


def matrix_stats(headers_list, distmats):
    """does, mean, median, stdev on a series of (dis)similarity matrices
    
    takes a series of parsed matrices (list of headers, list of numpy 2d arrays)
    headers must are either row or colunm headers (those must be identical)
    outputs headers (list), means, medians, stdevs (all numpy 2d arrays)
    """
    
    if len(set(map(tuple,headers_list))) > 1:
        print "error, not all input matrices have identical column/row headers"+\
          set(headers_list)
        sys.exit(1)
        
    all_mats = numpy.array(distmats) # 3d numpy array: mtx, row, col
    means = numpy.mean(all_mats, axis=0)
    medians = numpy.median(all_mats, axis=0)
    stdevs = numpy.std(all_mats, axis=0)
    
    return headers_list[0], means, medians, stdevs
  
  
usage_str = """usage: %prog [options] {-i INPUT_DIR -o OUTPUT_DIR}

[] indicates optional input (order unimportant) 
{} indicates required input (order unimportant) 

read in all (dis)similarity matrices in input_dir, write mean, median, stdev
to output folder.

input_dir must contain only (dis)similarity matrices, and only those you wish
to obtain stats on

Example usage:
python %prog -i dists/ -o dist_stats/

# Get detailed usage information
python pick_otus.py -h

"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    parser = OptionParser(usage=usage)

    parser.add_option('-i','--input_dir',\
          help='Path to input directory [REQUIRED]')

    parser.add_option('-o','--output_dir',\
          help='Path to store '+\
          'result files [REQUIRED]')

    opts,args = parser.parse_args()

    required_options = ['input_dir', 'output_dir']

    return opts,args
    
if __name__ == '__main__':
    options, args = parse_command_line_parameters()   
    indir = options.input_dir
    outdir = options.output_dir
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    #input    
    file_names = os.listdir(indir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    distmats = []
    headers_list = []
    for fname in file_names:
        f = open(os.path.join(indir,fname), 'U')
        headers, data = parse_distmat(f)
        f.close()
        distmats.append(data)
        headers_list.append(headers)
        
    #calcs
    headers, means, medians, stdevs = matrix_stats(headers_list, distmats)
    
    #output
    f = open(os.path.join(outdir,'means.txt'), 'w')
    f.write(format_distance_matrix(headers,means))
    f.close()

    f = open(os.path.join(outdir,'medians.txt'), 'w')
    f.write(format_distance_matrix(headers,medians))
    f.close()

    f = open(os.path.join(outdir,'stdevs.txt'), 'w')
    f.write(format_distance_matrix(headers,stdevs))
    f.close()