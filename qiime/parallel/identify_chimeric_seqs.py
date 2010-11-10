#!/usr/bin/env python
# Author: Jens Reeder (jens.reeder@gmail.com)
# identify_chimeric_seqs.py

from __future__ import division
from qiime.parallel.poller import basic_process_run_results_f
from qiime.parallel.util import get_rename_command

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso","Jens Reeder"] 
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

def get_job_commands(python_exe_fp, identify_chimeric_seqs_fp, fasta_fps,
                     output_dir, ref_seqs_fp, job_prefix, working_dir,
                     aligned_reference_seqs_fp, blast_db,
                     chimera_detection_method, min_div_ratio, num_fragments,
                     taxonomy_depth, max_e_value, id_to_taxonomy_fp,
                     command_prefix='', command_suffix=''):
#                     command_prefix='/bin/bash; ', command_suffix='; exit'):
    """Generate identify_chimeric_seqs commands which should be run
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_chimeric.txt']
    
    # Create lists to store the results
    commands = []
    result_filepaths = []
    
    # Iterate over the input files
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
        [fn % i for fn in out_filenames], working_dir, output_dir)
        result_filepaths += current_result_filepaths

        #Need to be filled
        optional_options = ""

        if chimera_detection_method=='blast_fragments':
            
            if ref_seqs_fp:
                optional_options += " -r %s" % ref_seqs_fp
            if blast_db:
                optional_options += " -b %s" % blast_db

            command = \
                '%s %s %s -i %s -t %s -m blast_fragments -o %s -n %s -d %s -e %s %s %s %s' %\
                (command_prefix,
                 python_exe_fp,
                 identify_chimeric_seqs_fp,
                 fasta_fp,
                 id_to_taxonomy_fp,
                 working_dir+"/"+out_filenames[0] % i,
                 num_fragments,
                 taxonomy_depth,
                 max_e_value,
                 optional_options,  
                 rename_command,
                 command_suffix)
            
        elif chimera_detection_method=='ChimeraSlayer':
            if min_div_ratio:
                optional_options += " --min_div_ratio %s" % min_div_ratio
            if ref_seqs_fp:
                optional_options += " -r %s" % ref_seqs_fp
            command = \
                '%s %s %s -i %s -a %s -m ChimeraSlayer -o %s %s %s %s' %\
                (command_prefix,
                 python_exe_fp,
                 identify_chimeric_seqs_fp,
                 fasta_fp,
                 aligned_reference_seqs_fp,
                 working_dir+"/"+out_filenames[0] % i,
                 optional_options,    
                 rename_command,
                 command_suffix)
        else:
           raise NotImplementedError
        commands.append(command)

    return commands, result_filepaths

def get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
    merge_map_filepath,deletion_list_filepath,process_run_results_f,\
    seconds_to_sleep,command_prefix='/bin/bash; ',command_suffix='; exit'):
    """Generate command to initiate a poller to monitior/process completed runs
    """
    
    result = '%s %s %s -f %s -p %s -m %s -d %s -t %d %s' % \
     (command_prefix,
      python_exe_fp,
      poller_fp,
      expected_files_filepath,
      process_run_results_f,
      merge_map_filepath,
      deletion_list_filepath,
      seconds_to_sleep,
      command_suffix)
      
    return result, []
