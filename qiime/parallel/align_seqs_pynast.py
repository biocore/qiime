#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# parallel_pynast.py

from __future__ import division
from qiime.parallel.util import get_rename_command

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

def get_job_commands(python_exe_fp,align_seqs_fp,fasta_fps,template_aln_fp,\
    pairwise_alignment_method,output_dir,blast_db,\
    min_length,min_percent_id,job_prefix,working_dir,command_prefix=None,\
    command_suffix=None):
    """Generate PyNAST commands which should be submitted to cluster
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_aligned.fasta', 
                     job_prefix + '.%d_failures.fasta',
                     job_prefix + '.%d_log.txt']
    
    # Initialize the command_prefix and command_suffix
    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    # Create lists to store the results
    commands = []
    result_filepaths = []
    
    # If there is a value for blast_db, pass it. If not, it
    # will be created on-the-fly. Note that on-the-fly blast dbs
    # are created with a string of random chars in the name, so this is safe.
    # They shouldn't overwrite one another, and will be cleaned up.
    if blast_db:
        blast_option = '-d %s' % blast_db
    else:
        blast_option = ''
    
    # Iterate over the input files
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [fn % i for fn in out_filenames],working_dir,output_dir)
        result_filepaths += current_result_filepaths
            
        command = \
         '%s %s %s %s -p %1.2f -e %d -m pynast -t %s -a %s -o %s -i %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          align_seqs_fp,\
          blast_option,\
          min_percent_id,\
          min_length,\
          template_aln_fp,\
          pairwise_alignment_method,
          working_dir,
          fasta_fp,
          rename_command,
          command_suffix)
          
        commands.append(command)

    return commands, result_filepaths
    