#!/usr/bin/env python
# File created on 20 Dec 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from optparse import make_option
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,make_option)
from qiime.filter import (filter_samples_from_distance_matrix, 
                         sample_ids_from_metadata_description,
                         get_seqs_to_keep_lookup_from_seq_id_file)
from qiime.parse import parse_otu_table, parse_distmat

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("",
 "Filter samples ids listed in sample_id_list.txt from dm.txt",
 "filter_distance_matrix.py -i dm.txt -o dm_out.txt -s sample_id_list.txt"),
 ("",
 "Filter samples ids in otu_table.txt from dm.txt",
 "filter_distance_matrix.py -i dm.txt -o dm_out.txt -t otu_table.txt")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_distance_matrix',
             help='the input distance matrix',type="existing_filepath"),
 make_option('-o','--output_distance_matrix',
             help='path to store the output distance matrix',type="new_filepath")]

script_info['optional_options'] = [
 make_option('--sample_id_fp',type="existing_filepath",
  help='A list of sample identifiers (or tab-delimited lines with'
  ' a sample identifier in the first field) which should be retained'),
 make_option('-t','--otu_table_fp',
  type="existing_filepath",help='the otu table filepath'),
 make_option('-m','--mapping_fp',help='path to the mapping file',type="existing_filepath"),\
 make_option('-s','--valid_states',help="string containing valid states, e.g. 'STUDY_NAME:DOG'"),
 make_option('--negate',default=False,
             action='store_true',
             help="discard specified samples (instead of keeping them) [default: %default]")]
script_info['version'] = __version__

def sample_ids_from_otu_table(otu_table_f):
    """ """
    return parse_otu_table(otu_table_f)[0]

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    output_f = open(opts.output_distance_matrix,'w')
    if opts.otu_table_fp:
        samples_to_keep = \
         sample_ids_from_otu_table(open(opts.otu_table_fp,'U'))
    elif opts.sample_id_fp:
        samples_to_keep = \
         get_seqs_to_keep_lookup_from_seq_id_file(open(opts.sample_id_fp,'U'))
    elif opts.mapping_fp and opts.valid_states:
        samples_to_keep = sample_ids_from_metadata_description(open(opts.mapping_fp,'U'),opts.valid_states)
    else:
        option_parser.error('must pass either --sample_id_fp, -t, or -m and -s')
    # note that negate gets a little weird here. The function we're calling removes the specified 
    # samples from the distance matrix, but the other QIIME filter scripts keep these samples specified. 
    # So, the interface of this script is designed to keep the specified samples, and therefore
    # negate=True is passed to filter_samples_from_distance_matrix by default.
    d = filter_samples_from_distance_matrix(
                               parse_distmat(open(opts.input_distance_matrix,'U')),
                               samples_to_keep,
                               negate=not opts.negate)
    output_f.write(d)
    output_f.close()
    


if __name__ == "__main__":
    main()