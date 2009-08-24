#!/usr/bin/env python
import numpy
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from pipe454.format import format_coords
from pipe454.parse import parse_distmat

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
    idxs_descending = eigvals.argsort()[::-1]
    coords = coords[idxs_descending]
    eigvals = eigvals[idxs_descending]
    pcnts = (numpy.abs(eigvals)/sum(numpy.abs(eigvals)))*100
    return format_coords(samples, coords.T, eigvals, pcnts)

def make_cmd_parser():
    """returns command-line options"""
    usage = 'usage: %prog dist_mtx_filepath [options]\n' +\
        'Default is to write to stdout'
    parser = OptionParser(usage=usage)
    parser.add_option('-o', '--out_path', dest='output_path', default=None,
        help='output path')   

    options, args = parser.parse_args()
    if (len(args) != 1 ):
        parser.error("incorrect number of arguments")
    return options, args

if __name__ == '__main__':
    options, args = make_cmd_parser()

    infilepath = args[0]
    f = open(infilepath)
    pcoa_res_string = pcoa(f)
    f.close()
    if options.output_path:
        f = open(options.output_path, 'w')
        f.write(pcoa_res_string)
        f.close()
    else:
        print pcoa_res_string
