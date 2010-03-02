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

script_info={}
script_info['brief_description']="""Summarize Taxa"""
script_info['script_description']="""The summarize_taxa.py script provides summary information of the representation of taxonomic groups within each sample. It takes an OTU table that contains taxonomic information as input. The taxonomic level for which the summary information is provided is designated with the -L option. The meaning of this level will depend on the format of the taxon strings that are returned from the taxonomy assignment step. The taxonomy strings that are most useful are those that standardize the taxonomic level with the depth in the taxonomic strings. For instance, for the RDP classifier taxonomy, Level 2 = Domain (e.g. Bacteria), 3 = Phylum (e.g. Firmicutes), 4 = Class (e.g. Clostridia), 5 = Order (e.g. Clostridiales), 6 = Family (e.g. Clostridiaceae), and 7 = Genus (e.g. Clostridium). By default, the relative abundance of each taxonomic group will be reported, but the raw counts can be returned if -r is set as False."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""The following command can be used to summarize taxa based on the Class, where the default parameters are used (no mapping file, delimiter for RDP ("-d ;") and output relative abundance ("-r True")) and the results are written to the file "Class.txt":""","""summarize_taxa.py -i otu_table.txt -L 4 -o Class.txt"""))
script_info['script_usage'].append(("""""","""Optionally the user can have the relative abundances added to the user-generated mapping file, by using the following command:""","""summarize_taxa.py -i otu_table.txt -L 4 -m Mapping_file.txt"""))
script_info['script_usage'].append(("""""","""Alternatively, the user may want to output the raw counts of each lineage within a sample, which can be used in the next step for making pie charts, by using the following command:""","""summarize_taxa.py -i otu_table.txt -L 4 -r False"""))
script_info['output_description']="""There are two possible output formats depending on whether or not a category mapping file is provided with the -m option. If a category mapping file is not provided, a table is returned where the taxonomic groups are each in a row and there is a column for each sample. If a category mapping file is provided, the summary information will be appended to this file. Specifically, a new column will be made for each taxonomic group to which the relative abundances or raw counts will be added to the existing rows for each sample. The addition of the taxonomic information to the category mapping file allows for taxonomic coloration of Principal coordinates plots in the 3d viewer. As described in the make_3d_plots.py section, principal coordinates plots can be dynamically colored based on any of the metadata columns in the category mapping file. Dynamic coloration of the plots by the relative abundances of each taxonomic group can help to distinguish which taxonomic groups are driving the clustering patterns.
"""

script_info['required_options']= [\
make_option('-i','--otu_file',action='store',\
          type='string',dest='otu_fp',help='Path to read '+\
          'otu file [REQUIRED]')
]

script_info['optional_options'] = [\
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

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

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
