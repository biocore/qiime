#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# rarefaction.py

from __future__ import division
from qiime.parallel.util import get_rename_command

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso","Justin Kuczynski"] 
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


def get_job_commands(python_exe_fp,rarefaction_fp,job_prefix,\
    input_fp,output_dir,working_dir,min_seqs,max_seqs,step,num_reps,\
    lineages_included, command_prefix=None,command_suffix=None):
    """Generate alpha diversity commands to be submitted to cluster
    """
    # Create data for each run (depth, output_fn)
    run_parameters = []
    for num_seqs in range(min_seqs,max_seqs+1, step):
        for rep_num in range(num_reps):
            run_parameters.append((\
             num_seqs,'rarefaction_%d_%d.txt' % (num_seqs,rep_num)))

    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    commands = []
    result_filepaths = []
    
    if lineages_included:
        lineages_included_param = '--lineages_included'
    else:
        lineages_included_param = ''
    
    for depth,output_fn in run_parameters:
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [output_fn],working_dir,output_dir)
        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -i %s -o %s %s -d %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          rarefaction_fp,\
          input_fp,
          working_dir + '/' + output_fn,
          lineages_included_param,
          depth,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths
