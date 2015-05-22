#!/usr/bin/env python
from __future__ import division
__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""takes a tree and bootstrap support file and writes a pdf, colored by
bootstrap support
"""
from matplotlib import use
use('Agg', warn=False)
from cogent.draw.dendrogram import SquareDendrogram
import os.path
import sys


def write_pdf_bootstrap_tree(tree, output_f, hits_dict):

    def f(node):
        if not node.Name:
            return 'black'
        tip_id = node.Name.split('/')[0]
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

    t = SquareDendrogram(tree)
    # Make output size proportional to the tree size.
    width = 8 * len(tree.tips())
    height = 8 * len(tree.tips())
    if width < 700:
        width = 700
    if height < 700:
        height = 700
    t.drawToPDF(output_f, width, height, edge_color_callback=f)
