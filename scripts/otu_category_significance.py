#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Doug Wendel"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Catherine Lozupone", "Jesse Stombaugh", "Doug Wendel"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Doug Wendel"
__email__ = "wendel@colorado.edu"
__status__ = "Pre-release"
 

from qiime.otu_category_significance import G_test_wrapper, ANOVA_wrapper
from qiime.util import parse_command_line_parameters
from optparse import make_option
from optparse import OptionParser
from os.path import split, splitext
from numpy import argsort
from cogent.util.dict2d import Dict2D
from cogent.maths.stats.test import calc_contingency_expected, G_fit_from_Dict2D,\
    ANOVA_one_way
from cogent.maths.stats.util import Numbers
from numpy import array
import sys

script_description = """Look for OTUs that are associated with a category."""

script_usage = """
Look for OTUs that are associated with a category. Currently can do:
1) perform g-test of independence to determine whether OTU presence
absence is associated with a category in the mapping file. Can
be used to look for co-occurrence using qPCR data or to look for 
associations with a categorical environmental variable.
2) perform ANOVA to determine whether OTU abundance is significantly
different across a category

python ~/repo/Qiime/qiime/OTU_category_significance.py -i otu_table.txt, -m category_mapping.txt -s g_test -f 10 -c category name -o output_fp -t None
"""

required_options = [\
    make_option('-i','--otu_table_fp', dest='otu_table_fp',\
        help='path to the otu table'),\
    make_option(    '-m','--category_mapping_fp',\
        dest='category_mapping_fp',\
        help='path to category mapping file'),\
    make_option('-c','--category', dest='category',\
        help='name of category over which to run the analysis')
]

optional_options = [\
    make_option('-s','--test', dest='test', default='g_test',\
        help='the type of statistical test to run. options are: ' +\
        'g_test: g test of independence: determines whether OTU ' +\
        'presence/absence is associated with a category ' +\
        'ANOVA: determines whether OTU abundance is associated with a ' +\
        'category'),\
    make_option('-o','--output_fp', dest='output_fp', \
        default= 'otu_category_G_test_results.txt',\
        help='path to output file. otu_category_G_test_results.txt by default'),\
    make_option('-f','--filter', dest='filter',\
        default= 10, \
        help='minimum number of samples that must contain the OTU for the ' +\
        'OTU to be included in the analysis. default value=10.'),\
    make_option('-t','--threshold', dest='threshold', default=None, \
        help='threshold under which to consider something absent: ' +\
        'Only used if you have numerical data that should be converted to ' +\
        'present or absent based on a threshold. Should be None for ' +\
        'categorical data. default value is None')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
        script_description=script_description,
        script_usage=script_usage,
        version=__version__,
        required_options=required_options,
        optional_options=optional_options)

    verbose = opts.verbose

    otu_table_fp = opts.otu_table_fp
    otu_table = open(otu_table_fp)
    output_fp = opts.output_fp

    category_mapping_fp = opts.category_mapping_fp
    category_mapping = open(category_mapping_fp,'U')

    filter = opts.filter
    category = opts.category
    threshold = opts.threshold
    if threshold and threshold != 'None':
        threshold = float(threshold)
    test = opts.test

    if test == 'g_test':
        G_test_wrapper(otu_table, category_mapping, category, threshold, \
            filter, output_fp)
    elif test == 'ANOVA':
        ANOVA_wrapper(otu_table, category_mapping, category, threshold, \
            filter, output_fp)
      

if __name__ == "__main__":
    main()