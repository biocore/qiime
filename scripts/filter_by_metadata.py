#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"	
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Antonio Gonzalez Pena"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from string import strip
from qiime.filter_by_metadata import filter_otus_and_map

script_description = """This filter allows for the removal of sequences and OTUs that either do or don't match specified 
metadata, for instance, isolating samples from a specific set of studies or body sites. This script 
identifies samples matching the specified metadata criteria, and outputs a filtered mapping file 
and OTU table containing only the specified samples."""

script_usage = """ The following command can be used, where all options are passed (using the resulting OTU file from 
  pick_otus.py, the original Fasting_Map file, and keeping only the Control sequences in the 
  Treatment field) with the resulting data being written to seqs_otus.txt.filtered.xls and 
  Fasting_Map.txt.filtered.xls:

  filter_otus_by_sample.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control' 

  Some variations (not so useful on this dataset, but more useful on larger datasets) are:
  - Keeping both Control and Fast in the Treatment field (i.e. keeping everything): 
       filter_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control,Fast'
  - Excluding Fast in the Treatment field (same as the first example) - the syntax here is * to keep 
    everything, then !Fast to eliminate the Fast group:
       filter_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:*,!Fast'
  - Keeping only samples with both Control in the Treatment field and 20061218 in the DOB field:
       filter_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control,DOB: 20061218'
  - Keeping only samples with Control in the Treatment field and OTUs with counts of at least 
    5 across samples: 
       filter_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control' -n 5
       
  Note that the filtered mapping file will automatically exclude any columns that are the same for 
  all the samples that are left, and will also exclude (except for SampleID) any columns that are 
  different for all the samples that are left, making it more useful for downstream analyses with the 
  coloring tools.
"""


required_options = [\
 make_option('-i', '--otu', dest='otu_fname',\
        help='name of otu file [REQUIRED]'),\
 make_option('-m', '--map', dest='map_fname',\
        help='name of map file [REQUIRED]'),\
 make_option('-s', '--states', dest='valid_states',\
        help="string containing valid states, e.g. 'STUDY_NAME:DOG'")
]

optional_options = [\
 make_option('-o', '--otu_outfile', dest='otu_out_fname', default=None,\
        help='name of otu output file, default is otu_filename.filtered.xls'),\
 make_option('-p', '--map_outfile', dest='map_out_fname', default=None,\
        help='name of map output file, default is map_filename.filtered.xls'),\
 make_option('-n', '--num_seqs_per_otu', dest='num_seqs_per_otu',\
        type=int, default=1, help='minimum counts across samples to keep OTU,' +\
        ' default is only to keep OTUs that are present in the samples.')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    map_file_name, otu_file_name, valid_states_str, num_seqs_per_otu = \
        opts.map_fname, opts.otu_fname, opts.valid_states, \
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