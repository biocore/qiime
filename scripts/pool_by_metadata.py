#!/usr/bin/env python
# File created on 18 Oct 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from string import strip
from qiime.filter_by_metadata import parse_metadata_state_descriptions,\
get_sample_ids
from qiime.pool_by_metadata import pool_map, pool_otu_table
from qiime.parse import parse_mapping_file, parse_otu_table

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""pool samples in OTU table and mapping file based on sample metadata from mapping file"""
script_info['script_description']="""this script outputs a new otu table and mapping file with some samples removed and replaced with one pooled sample. the new pooled sample will have fields in the mapping file the same as its constituent samples, if all are idential. Else it will just say 'multipleValues'."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""The following command pools all the Control samples into one sample named 'pooledControl'. The resulting data is written to seqs_otus.txt.filtered.xls and Fasting_Map.txt.filtered.xls:""","""%prog -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control' -l pooledControl"""))
script_info['script_usage'].append(("""""","""Some variations are: 

Pooling all samples in both Control and Fast in the Treatment field (i.e. pooling everything):""","""%prog -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control,Fast' -l pooledSamples"""))
script_info['script_usage'].append(("""""","""Excluding Fast in the Treatment field - the syntax here is "*" to keep everything, then !Fast to eliminate the Fast group:""","""%prog -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:*,!Fast' -l pooledNonFast"""))


script_info['output_description']="""The result is a pooled OTU table and mapping file."""
script_info['required_options']=[\
   options_lookup['otu_table_as_primary_input'],\
 make_option('-m', '--map', dest='map_fname',\
        help='path to the map file [REQUIRED]'),\
 make_option('-s', '--states', dest='valid_states',\
        help="string containing valid states, e.g. 'STUDY_NAME:DOG'")\
]
script_info['optional_options']=[\
 make_option('-o', '--otu_outfile', dest='otu_out_fname', default=None,\
        help='name of otu output file, default is otu_filename.pooled.txt'),\
 make_option('-p', '--map_outfile', dest='map_out_fname', default=None,\
        help='name of map output file, default is map_filename.pooled.txt'),\
 make_option('-l', '--pooled_sample_name', dest='pooled_sample_name', default='pooledSample',\
        help='new sample name used in new mapping file and new otu table'),\

]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    map_file_name, otu_file_name, valid_states_str = \
        opts.map_fname, opts.otu_table_fp, opts.valid_states
    map_infile = open(map_file_name, 'U')
    otu_infile = open(otu_file_name, 'U')

    map_out_fname = opts.map_out_fname
    otu_out_fname = opts.otu_out_fname
    
    if map_out_fname is None:
        map_out_fname = map_file_name + '.pooled.txt'

    if otu_out_fname is None:
        otu_out_fname = otu_file_name + '.pooled.txt'

    # write out the filtered mapping file
    map_outfile = open(map_out_fname, 'w')
    otu_outfile = open(otu_out_fname, 'w')

    map_data, map_header, map_comments = parse_mapping_file(map_infile)
    map_infile.close()
    map_infile = open(map_file_name, 'U') # reopen for later
    valid_states = parse_metadata_state_descriptions(valid_states_str)
    sample_ids_to_pool = get_sample_ids(map_data, map_header, valid_states)
    
    pool_map(map_infile, map_outfile,
        opts.pooled_sample_name, sample_ids_to_pool)
    pool_otu_table(otu_infile, otu_outfile,
        opts.pooled_sample_name, sample_ids_to_pool)

if __name__ == "__main__":
    main()