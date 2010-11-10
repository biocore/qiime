#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski",\
        "Julia Goodrich", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, \
        convert_otu_table_relative
from optparse import make_option
from qiime.summarize_taxa import make_summary, add_summary_mapping
from sys import stdout, stderr
from qiime.parse import parse_otu_table, parse_mapping_file
from qiime.format import write_summarize_taxa, write_add_taxa_summary_mapping,\
        format_summarize_taxa, format_add_taxa_summary_mapping

script_info={}
script_info['brief_description']="""Summarize Taxa"""
script_info['script_description']="""The summarize_taxa.py script provides summary information of the representation of taxonomic groups within each sample. It takes an OTU table that contains taxonomic information as input. The taxonomic level for which the summary information is provided is designated with the -L option. The meaning of this level will depend on the format of the taxon strings that are returned from the taxonomy assignment step. The taxonomy strings that are most useful are those that standardize the taxonomic level with the depth in the taxonomic strings. For instance, for the RDP classifier taxonomy, Level 2 = Domain (e.g. Bacteria), 3 = Phylum (e.g. Firmicutes), 4 = Class (e.g. Clostridia), 5 = Order (e.g. Clostridiales), 6 = Family (e.g. Clostridiaceae), and 7 = Genus (e.g. Clostridium). By default, the relative abundance of each taxonomic group will be reported, but the raw counts can be returned if -a is passed."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""The following command can be used to summarize taxa based on the Class, where the default parameters are used (no mapping file, delimiter for RDP ("-d ;") and output relative abundance) and the results are written to the file "Class.txt":""","""summarize_taxa.py -i otu_table.txt -L 4 -o Class.txt"""))
script_info['script_usage'].append(("""""","""Optionally the user can have the relative abundances added to the user-generated mapping file, by using the following command:""","""summarize_taxa.py -i otu_table.txt -L 4 -m Mapping_file.txt"""))
script_info['script_usage'].append(("""""","""Alternatively, the user may want to output the raw counts of each lineage within a sample, which can be used in the next step for making pie charts, by using the following command:""","""summarize_taxa.py -i otu_table.txt -L 4 -a"""))
script_info['output_description']="""There are two possible output formats depending on whether or not a mapping file is provided with the -m option. If a mapping file is not provided, a table is returned where the taxonomic groups are each in a row and there is a column for each sample. If a mapping file is provided, the summary information will be appended to this file. Specifically, a new column will be made for each taxonomic group to which the relative abundances or raw counts will be added to the existing rows for each sample. The addition of the taxonomic information to the mapping file allows for taxonomic coloration of Principal coordinates plots in the 3d viewer. As described in the make_3d_plots.py section, principal coordinates plots can be dynamically colored based on any of the metadata columns in the mapping file. Dynamic coloration of the plots by the relative abundances of each taxonomic group can help to distinguish which taxonomic groups are driving the clustering patterns.
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
make_option('-m','--mapping',action='store',\
          type='string',dest='mapping', 
          help='if supplied - the taxon information will be added ' +\
             'to the mapping file. This mapping file can ' +\
             'be used to color PCoA plots by taxon abundance or to ' +\
             'perform statistical tests of taxon/mappingy associations.'),
make_option('-d','--delimiter',action='store',\
          type='string',dest='delimiter',default=';', 
          help='Delimitor that separates taxonomy categories.[default: %default]'),
make_option('-r', '--relative_abundance', action='store',\
        dest='relative_abundance', default='', \
        help='DEPRECATED: please use -a/--absolute_abundance to disable ' +\
        'relative abundance [default: %default]'),
make_option('-a', '--absolute_abundance', action='store_true',\
        dest='absolute_abundance', default=False, \
        help='If present, reports the absolute abundance of the lineage in ' +\
            'each sample. By default uses relative abundance [default: %default]')
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_fname = opts.out_fp
    otu_fp = opts.otu_fp
    otu_table = parse_otu_table(open(otu_fp, 'U'))
    delimiter = opts.delimiter
    mapping = opts.mapping
    level = opts.level

    if mapping:
        mapping_file = open(mapping, 'U')
        mapping, header, comments = parse_mapping_file(mapping_file)

    if opts.relative_abundance != '':
        raise option_parser.error("Depreciated. Please use --absolute_abundances to disable relative abundance")

    if not opts.absolute_abundance:
        otu_table = convert_otu_table_relative(otu_table)

    if output_fname:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout

    if mapping:
        summary, tax_order = add_summary_mapping(otu_table, mapping, level)
        if output_fname:
            write_add_taxa_summary_mapping(summary,tax_order,mapping,header,output_fname,delimiter)
        else:
            print ''.join(format_add_taxa_summary_mapping(summary,tax_order,mapping,header,delimiter))

    else:
        summary, header = make_summary(otu_table, level)
        if output_fname:
            write_summarize_taxa(summary, header, output_fname, delimiter)
        else:
            print ''.join(format_summarize_taxa(summary,header,delimiter))


if __name__ == "__main__":
    main()
