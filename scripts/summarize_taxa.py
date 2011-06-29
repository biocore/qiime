#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski",\
        "Julia Goodrich", "Daniel McDonald", "Antonio Gonzalez Pena",
        "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters, \
        convert_otu_table_relative
from qiime.util import make_option,get_options_lookup,create_dir
from qiime.summarize_taxa import make_summary, add_summary_mapping
from sys import stdout, stderr
from qiime.parse import parse_otu_table, parse_mapping_file
from qiime.format import write_summarize_taxa, write_add_taxa_summary_mapping,\
        format_summarize_taxa, format_add_taxa_summary_mapping
from os.path import split,splitext,join

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Summarize Taxa"""
script_info['script_description']="""The summarize_taxa.py script provides summary information of the representation of taxonomic groups within each sample. It takes an OTU table that contains taxonomic information as input. The taxonomic level for which the summary information is provided is designated with the -L option. The meaning of this level will depend on the format of the taxon strings that are returned from the taxonomy assignment step. The taxonomy strings that are most useful are those that standardize the taxonomic level with the depth in the taxonomic strings. For instance, for the RDP classifier taxonomy, Level 2 = Domain (e.g. Bacteria), 3 = Phylum (e.g. Firmicutes), 4 = Class (e.g. Clostridia), 5 = Order (e.g. Clostridiales), 6 = Family (e.g. Clostridiaceae), and 7 = Genus (e.g. Clostridium). By default, the relative abundance of each taxonomic group will be reported, but the raw counts can be returned if -a is passed."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""The following command can be used to summarize taxa based on the Class, where the default parameters are used (no mapping file, delimiter for RDP ("-d ;") and output relative abundance) and the results are written to the directory "Class":""","""summarize_taxa.py -i otu_table.txt -L 4 -o ./Class"""))
script_info['script_usage'].append(("""""","""Optionally the user can have the relative abundances added to the user-generated mapping file, by using the following command:""","""summarize_taxa.py -i otu_table.txt -L 4 -m Mapping_file.txt"""))
script_info['script_usage'].append(("""""","""Alternatively, the user may want to output the raw counts of each lineage within a sample, which can be used in the next step for making area, bar and pie charts, by using the following command:""","""summarize_taxa.py -i otu_table.txt -L 4 -a"""))
script_info['output_description']="""There are two possible output formats depending on whether or not a mapping file is provided with the -m option. If a mapping file is not provided, a table is returned where the taxonomic groups are each in a row and there is a column for each sample. If a mapping file is provided, the summary information will be appended to this file. Specifically, a new column will be made for each taxonomic group to which the relative abundances or raw counts will be added to the existing rows for each sample. The addition of the taxonomic information to the mapping file allows for taxonomic coloration of Principal coordinates plots in the 3d viewer. As described in the make_3d_plots.py section, principal coordinates plots can be dynamically colored based on any of the metadata columns in the mapping file. Dynamic coloration of the plots by the relative abundances of each taxonomic group can help to distinguish which taxonomic groups are driving the clustering patterns.
"""

script_info['required_options']= [\
    make_option('-i','--otu_table_fp', dest='otu_table_fp',
        help='Input OTU table filepath [REQUIRED]',
        type='existing_filepath'),
]
script_info['optional_options'] = [\
    make_option('-L','--level',default='2,3,4,5,6' ,
        help='Taxonomic level to summarize by. [default: %default]'),
    make_option('-m','--mapping', 
        help='Input metadata mapping filepath. If supplied, then the taxon' +\
        ' information will be added to this file. This option is ' +\
        ' useful for coloring PCoA plots by taxon abundance or to ' +\
        ' perform statistical tests of taxon/mapping associations.',
        type='existing_filepath'),
    make_option('-d','--delimiter',action='store',type='string',
        dest='delimiter',default=';', 
        help='Delimitor separating taxonomy levels. [default: %default]'),
    make_option('-r', '--relative_abundance', action='store',\
        dest='relative_abundance', default='', \
        help='DEPRECATED: please use -a/--absolute_abundance to disable ' +\
        'relative abundance [default: %default]'),
    make_option('-a', '--absolute_abundance', action='store_true',\
        dest='absolute_abundance', default=False, \
        help='If present, the absolute abundance of the lineage in ' +\
        ' each sample is reported. By default, this script uses relative' +\
        ' abundance [default: %default]'),
    make_option('-l', '--lower_percentage', type='float', default=None, \
        help='If present, OTUs having higher absolute abundance are' +\
        ' trimmed. To remove OTUs that make up more than 5% of the total' +\
        ' dataset you would pass 0.05. [default: %default]'),
    make_option('-u', '--upper_percentage', type='float', default=None, \
        help='If present, OTUs having lower absolute abundance are' +\
        ' trimmed. To remove the OTUs that makes up less than 45% of the' +\
        ' total dataset you would pass 0.45. [default: %default]'),
    options_lookup['output_dir'],
]
script_info['option_label']={'otu_table_fp':'OTU table filepath',
                             'output_fp': 'Output filepath',
                             'mapping':'QIIME-formatted mapping filepath',
                             'level':'Summarize level',
                             'delimiter': 'Taxonomic delimiter',
                             'relative_abundance':'Use relative abundance',
                             'absolute_abundance':'Use absolute abundance',
                             'lower_percentage':'Top % of OTUs to remove',
                             'upper_percentage':'Bottom % of OTUs to remove'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    lower_percentage = opts.lower_percentage
    upper_percentage = opts.upper_percentage
    otu_table_fp = opts.otu_table_fp
    otu_table = parse_otu_table(open(otu_table_fp, 'U'))
    delimiter = opts.delimiter
    mapping_fp = opts.mapping
    levels = opts.level.split(',')

    if upper_percentage!=None and lower_percentage!=None:
        raise ValueError("upper_percentage and lower_percentage are mutually exclusive")
    
    if upper_percentage!=None and lower_percentage!=None and mapping:
        raise ValueError("upper_percentage and lower_percentage can not be using with mapping file")
        
    if upper_percentage!=None and (upper_percentage<0 or upper_percentage>1.0):
        raise ValueError('max_otu_percentage should be between 0.0 and 1.0')
    
    if lower_percentage!=None and (lower_percentage<0 or lower_percentage>1.0):
        raise ValueError('lower_percentage should be between 0.0 and 1.0')
        
    if mapping_fp:
        mapping_file = open(mapping_fp, 'U')
        mapping, header, comments = parse_mapping_file(mapping_file)
        
        # use the input Mapping file for producing the output filenames
        map_dir_path,map_fname=split(mapping_fp)
        map_basename,map_fname_ext=splitext(map_fname)

    if opts.relative_abundance != '':
        raise option_parser.error("Deprecated. Please use --absolute_abundances to disable relative abundance")

    if not opts.absolute_abundance:
        otu_table = convert_otu_table_relative(otu_table)

    # introduced output directory to will allow for multiple outputs
    if opts.output_dir:
        create_dir(opts.output_dir,False)
        output_dir_path=opts.output_dir
    else:
        output_dir_path='./'

    # use the input OTU table to produce the output filenames
    dir_path,fname=split(otu_table_fp)
    basename,fname_ext=splitext(fname)
    
    # Iterate over the levels and generate a summarized taxonomy for each
    for level in levels:
        if mapping_fp:
            #define output filename
            output_fname = join(output_dir_path,
                                        map_basename+'_L%s.txt' % (level))
                                        
            summary, tax_order = add_summary_mapping(otu_table, mapping,
                                                     int(level))
            write_add_taxa_summary_mapping(summary,tax_order,mapping,
                                            header,output_fname,delimiter)
        else:
            #define output filename
            output_fname = join(output_dir_path,basename+'_L%s.txt' % (level))
            
            summary, header = make_summary(otu_table, int(level),
                                            upper_percentage, lower_percentage)
            write_summarize_taxa(summary, header, output_fname, delimiter)
            

if __name__ == "__main__":
    main()
