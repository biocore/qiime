#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"	
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Antonio Gonzalez Pena","Greg Caporaso"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Tony Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from string import strip
from qiime.filter_by_metadata import filter_otus_and_map

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Filter OTU table by removal of specified metadata"""
script_info['script_description']="""This filter allows for the removal of sequences and OTUs that either do or don't match specified metadata, for instance, isolating samples from a specific set of studies or body sites. This script identifies samples matching the specified metadata criteria, and outputs a filtered mapping file and OTU table containing only the specified samples."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""The following command can be used, where all options are passed (using the resulting OTU file from make_otu_table.py, the original Fasting_Map.txt, and keeping only the Control sequences in the Treatment field) with the resulting data being written to otu_table.txt.filtered.xls and Fasting_Map.txt.filtered.xls:""","""%prog -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control'"""))
script_info['script_usage'].append(("""""","""Some variations (not so useful on this dataset, but more useful on larger datasets) are: 

Keeping both Control and Fast in the Treatment field (i.e. keeping everything):""","""%prog -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control,Fast'"""))
script_info['script_usage'].append(("""""","""Excluding Fast in the Treatment field (same as the first example) - the syntax here is "*" to keep everything, then !Fast to eliminate the Fast group:""","""%prog -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:*,!Fast'"""))
script_info['script_usage'].append(("""""","""Keeping only samples with both Control in the Treatment field and 20061218 in the DOB field:""","""        %prog -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control;DOB:20061218'"""))
script_info['script_usage'].append(("""""","""Keeping only samples with Control in the Treatment field and OTUs with counts of at least 5 across samples:""","""%prog -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control' -n 5"""))
script_info['script_usage'].append(("""""","""Note that the filtered mapping file will automatically exclude any columns that are the same for all the samples that are left, and will also exclude (except for SampleID) any columns that are different for all the samples that are left, making it more useful for downstream analyses with the coloring tools.""",""""""))
script_info['output_description']="""The result is a filtered OTU table and mapping file meeting the desired criteria."""
script_info['required_options']=[\
   options_lookup['otu_table_as_primary_input'],\
 make_option('-m', '--map', dest='map_fname',\
        help='path to the map file [REQUIRED]'),\
 make_option('-s', '--states', dest='valid_states',\
        help="string containing valid states, e.g. 'STUDY_NAME:DOG'")\
]
script_info['optional_options']=[\
 make_option('-o', '--otu_outfile', dest='otu_out_fname', default=None,\
        help='name of otu output file, default is otu_filename.filtered.xls'),\
 make_option('-p', '--map_outfile', dest='map_out_fname', default=None,\
        help='name of map output file, default is map_filename.filtered.xls'),\
 make_option('-n', '--num_seqs_per_otu', dest='num_seqs_per_otu',\
        type=int, default=1, help='minimum counts across samples to keep OTU,' +\
        ' default is only to keep OTUs that are present in the samples.')\
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    map_file_name, otu_file_name, valid_states_str, num_seqs_per_otu = \
        opts.map_fname, opts.otu_table_fp, opts.valid_states, \
        opts.num_seqs_per_otu

    map_infile = open(map_file_name, 'U')
    otu_infile = open(otu_file_name, 'U')

    map_out_fname = opts.map_out_fname
    otu_out_fname = opts.otu_out_fname
    
    if map_out_fname is None and otu_out_fname is None:
        map_out_fname = map_file_name + '.filtered.xls'
        otu_out_fname = otu_file_name + '.filtered.xls'
    elif otu_out_fname is None:
        otu_out_fname = map_out_fname + '-otu.filtered.xls'
    elif map_out_fname is None:
        map_out_fname = otu_out_fname + '-map.filtered.xls'

    # write out the filtered mapping file
    map_outfile = open(map_out_fname, 'w')
    otu_outfile = open(otu_out_fname, 'w')

    filter_otus_and_map(map_infile, otu_infile, map_outfile, otu_outfile,
        valid_states_str, num_seqs_per_otu)


if __name__ == "__main__":
    main()