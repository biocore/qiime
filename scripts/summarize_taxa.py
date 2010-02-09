#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.summarize_taxa import make_new_summary_file,add_summary_category_mapping
from sys import stdout, stderr

script_description = """Creates OTU tables that summarize table with taxa in \
last field."""

script_usage = """Creates an OTU table by taxa by a specified level. It uses \
the OTU file otus.txt (-i). The user can specify a taxonomic level (-L) the \
default is 2. They can also specify the delimitor for the taxonomy categories \
(-d) the default is ';' The output is an OTU table where each taxon \
(columns) was present in each sample (rows).  Samples present in category \
mapping file but absent from otu table are not included in output

python ~/code/Qiime/trunk/qiime/summarize_taxa.py -i otus.txt \
-o /Users/bob/qiime_run/
"""

required_options = [\
make_option('-i','--otu_file',action='store',\
          type='string',dest='otu_fp',help='Path to read '+\
          'otu file [REQUIRED]')
]

optional_options = [\
make_option('-o','--output_file',action='store',\
          type='string',dest='out_fp',help='Path to write '+\
          'output file'),
make_option('-L','--level',action='store',\
          type='int',dest='level', default=2, 
          help='Level of taxonomy to use [default: %default]'),
make_option('-m','--category_mapping',action='store',\
          type='string',dest='category_mapping', 
          help='if supplied - the taxon information will be added ' +\
             'to the category mapping file. This mapping file can ' +\
             'be used to color PCoA plots by taxon abundance or to ' +\
             'perform statistical tests of taxon/category associations.'),
make_option('-d','--delimitor',action='store',\
          type='string',dest='delimitor',default=';', 
          help='Delimitor that separates taxonomy categories.[default: %default]'),
make_option('-r', '--relative_abundance', action='store',\
        type='string', dest='relative_abundance', default='True', \
        help='If True, reports the relative abundance of the lineage in ' +\
            'each sample. If False, reports the raw counts. [default: %default]')
]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    output_fname = opts.out_fp
    otu_fp = opts.otu_fp
    otu_table = open(otu_fp, 'U')
    delimitor = opts.delimitor
    category_mapping = opts.category_mapping
    if category_mapping:
        category_mapping = open(category_mapping, 'U')
    level = opts.level
    relative_abundance=opts.relative_abundance
    if output_fname:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout
    if not category_mapping:
        output = make_new_summary_file(otu_table, level, delimitor, \
            relative_abundance)
    else:
        output = add_summary_category_mapping(otu_table, category_mapping, \
            level, delimitor, relative_abundance)
    outfile.write(''.join(output))


if __name__ == "__main__":
    main()
