#!/usr/bin/env python

from __future__ import division
from qiime.parallel.util import get_rename_command

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

def get_commands(python_exe_fp,assign_taxonomy_fp,confidence,job_prefix,\
    fasta_fps,rdp_jar_fp,output_dir,working_dir,\
    command_prefix=None,command_suffix=None,\
    id_to_taxonomy_fp=None,reference_seqs_fp=None):
    """Generate RDP classifier commands which should be submitted to cluster
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_tax_assignments.log', 
                     job_prefix + '.%d_tax_assignments.txt']
    
    command_prefix = command_prefix or\
     '/bin/bash; export RDP_JAR_PATH=%s; ' % rdp_jar_fp
    command_suffix = command_suffix or\
     '; exit'
    
    rdp_extra_params = ''
    if id_to_taxonomy_fp and reference_seqs_fp:
        rdp_extra_params = '-t %s -r %s' % (id_to_taxonomy_fp, reference_seqs_fp)
    
    commands = []
    result_filepaths = []
    
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [fn % i for fn in out_filenames],working_dir,output_dir)#,\
         #id_to_taxonomy_fp,reference_seqs_fp)
        result_filepaths += current_result_filepaths
        command = '%s %s %s %s -c %1.2f -m rdp -o %s -i %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          assign_taxonomy_fp,\
          rdp_extra_params,
          confidence,
          working_dir,
          fasta_fp,
          rename_command,
          command_suffix)
        commands.append(command)
        
    return commands, result_filepaths
    