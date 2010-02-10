#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters, matrix_stats
from optparse import make_option
import os
from qiime.parse import parse_distmat
from qiime.format import format_distance_matrix


script_description = """Description:
read in all (dis)similarity matrices in input_dir, write mean, median, stdev
to output folder.  outputs are in distance matrix format, each value is the 
(mean, median, or stdev) of that element in all the input distmats

input_dir must contain only (dis)similarity matrices, and only those you wish
to obtain stats on"""

script_usage = """python %prog -i dists/ -o dist_stats/"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
 make_option('-i','--input_dir',
       help='Path to input directory'),

 make_option('-o','--output_dir',
       help='Path to store result files')
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    
    indir = opts.input_dir
    outdir = opts.output_dir
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

if __name__ == "__main__":
    main()