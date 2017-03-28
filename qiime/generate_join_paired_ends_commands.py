#!/usr/bin/env python
from __future__ import division
 
__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import join, basename
                                
def create_commands(pairs, base_output_dir, optional_params = "",
        leading_text = "", trailing_text = "", include_input_dir_path=False,
        remove_filepath_in_name=False):
    """ Creates, calls commands for join_paired_ends.py

    pairs: dictionary of forward:reverse read filepaths
    base_output_dir: output directory to write log, stitched reads
    optional_params: added parameters to join_paired_ends.py calls
    leading_text: Text to add before join_paired_ends.py call
    trailing_text: Text to add after join_paired_ends.py call
    include_input_dir_path: If True, include input directory in output
        directory names
    remove_filepath_in_name: If True, the base filename will not be used in the
        output directory names.
    """
    
    commands = []
    for curr_fp in pairs:
        if include_input_dir_path:
            added_output_str = curr_fp.split(".fastq")[0].split('/')[-2]
        else:
            added_output_str = ""
        if not remove_filepath_in_name:
            added_output_str += basename(curr_fp).split('.fastq')[0]
            
        curr_outputdir = join(base_output_dir, added_output_str) 
        command = "%sjoin_paired_ends.py %s -f %s -r %s -o %s %s" %\
            (leading_text, optional_params, curr_fp, pairs[curr_fp],
            curr_outputdir, trailing_text)
            
        commands.append([('', command)])
                    
    return commands
    
def get_pairs(all_files, read1_indicator, read2_indicator):
    """ Finds pairs of files from a list of files
    
    all_files: list of filepaths
    read1_indicator: string indicating read 1 of a pair
    read2_indicator: string indicating read 2 of a pair
    """
    
    pairs = {}
    
    read1_files = []
    read2_files = []
    
    for curr_file in all_files:
        curr_file_string_r1 = curr_file.split(read1_indicator)
        curr_file_string_r2 = curr_file.split(read2_indicator)
        if len(curr_file_string_r1) == 2:
            read1_files.append(curr_file_string_r1)
        elif len(curr_file_string_r2) == 2:
            read2_files.append(curr_file_string_r2)
        else:
            raise ValueError,("Invalid filename found for splitting on input "+\
                "for file %s, " % curr_file + "check input read1_indicator "+\
                "and read2_indicator parameters as well.")
                
    for curr_read1 in read1_files:
        for curr_read2 in read2_files:
            if curr_read1 == curr_read2:
                pairs[read1_indicator.join(curr_read1)] =\
                    read2_indicator.join(curr_read2)
    
    return pairs
