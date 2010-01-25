#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
from cogent.draw.dendrogram import SquareDendrogram
from cogent import LoadTree
from qiime.parse import parse_bootstrap_support
import os.path
import sys

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"

"""takes a tree and bootstrap support file and writes a pdf, colored by
bootstrap support
"""


def main(opts):
    tree = LoadTree(opts.master_tree)
    support_file = open(opts.support)
    bootstraps = parse_bootstrap_support(support_file)
    support_file.close()
    write_pdf_bootstrap_tree(tree, opts.output_file, bootstraps)

            
usage_str = """usage: %prog [options] {-m MASTER_TREE -s SUPPORT -o OUTPUT_FILE}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
python %prog -m master_tree.tre -s jackknife_support.txt -o jackknife_samples.pdf

"""
def parse_command_line_parameters():
    """returns command-line options"""

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-m', '--master_tree',
        help='master tree filepath [REQUIRED]')

    parser.add_option('-s', '--support',
        help='path to bootstrap support file [REQUIRED]')

    parser.add_option('-o', '--output_file',
        help="output file name, will attach .pdf if not present [REQUIRED]")

    opts, args = parser.parse_args()
    if len(args) != 0:
        parser.error("positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')
         
    required_options = ['master_tree','support','output_file']
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
    return opts, args



def write_pdf_bootstrap_tree(tree, output_f, hits_dict):

    def f(node):
        if not node.Name:
            return 'black'
        tip_id=node.Name.split('/')[0]
        try:
            if hits_dict[tip_id] < .25:
                return 'blue'
            elif hits_dict[tip_id] < .5:
                return 'green'
            elif hits_dict[tip_id] < .75:
                return 'yellow'
            elif hits_dict[tip_id] <= 1.1:
                return 'red'
            return 'black'
        except:
            return 'black'
            
    t=SquareDendrogram(tree)
    # Make output size proportional to the tree size.
    width=8*len(tree.tips())
    height=8*len(tree.tips())
    if width<700:
        width=700
    if height<700:
        height=700
    if not output_f[-4:] == ".pdf":
        output_f = output_f + ".pdf"
    t.drawToPDF(output_f, width, height, edge_color_callback=f) 


if __name__ == '__main__':
    opts,args = parse_command_line_parameters()
    main(opts)
