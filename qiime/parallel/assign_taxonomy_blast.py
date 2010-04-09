#!/usr/bin/env python

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

def get_job_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,\
    e_value,blast_db,job_prefix,\
    blastmat_path,fasta_fps,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    """Generate BlastTaxonAssiger classifier commands to be submitted to cluster
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_tax_assignments.log', 
                     job_prefix + '.%d_tax_assignments.txt']

    command_prefix = command_prefix or\
     '/bin/bash; cd %s; export BLASTMAT=%s;' \
       % (working_dir,blastmat_path)
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    result_filepaths = []
    
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [fn % i for fn in out_filenames],working_dir,output_dir)
        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -o %s -m blast -e %s -b %s -i %s -t %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          assign_taxonomy_fp,\
          working_dir,
          e_value,
          blast_db,
          fasta_fp,
          id_to_taxonomy_fp,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths
