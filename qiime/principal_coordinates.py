#!/usr/bin/env python
import numpy
import sys
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from qiime.format import format_coords
from qiime.parse import parse_distmat

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"


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

usage_str = """usage: %prog [options] {-i INPUT_PATH -o OUTPUT_PATH}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

python %prog -i unifrac_dist_mtx.txt -o unifrac_pcoa.txt
"""
def parse_command_line_parameters():
    """returns command-line options"""
    
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
    
    parser.add_option('-i', '--input_path',
        help='input filepath. [REQUIRED]')
        
    parser.add_option('-o', '--output_path',
        help='output filepath. [REQUIRED]')

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
    opts, args = parse_command_line_parameters()

    infilepath = opts.input_path
    f = open(infilepath,'U')
    pcoa_res_string = pcoa(f)
    f.close()

    f = open(opts.output_path, 'w')
    f.write(pcoa_res_string)
    f.close()
