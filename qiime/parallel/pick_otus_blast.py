#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# pick_otus_blast.py

from __future__ import division
from qiime.parallel.poller import basic_process_run_results_f
from qiime.parallel.util import get_rename_command

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

def get_job_commands(python_exe_fp,pick_otus_fp,fasta_fps,\
    output_dir,blast_db,job_prefix,working_dir,max_e_value,\
    similarity,command_prefix='/bin/bash; ',command_suffix='; exit'):
    """Generate pick_otus commands which should be submitted to cluster
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_otus.log', 
                     job_prefix + '.%d_otus.txt']
    
    # Create lists to store the results
    commands = []
    result_filepaths = []
    
    # Iterate over the input files
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [fn % i for fn in out_filenames],working_dir,output_dir)
        result_filepaths += current_result_filepaths
            
        command = \
         '%s %s %s -i %s -b %s -m blast -o %s -e %s -s %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          pick_otus_fp,\
          fasta_fp,\
          blast_db,\
          working_dir,\
          max_e_value,\
          similarity,\
          rename_command,
          command_suffix)
          
        commands.append(command)

    return commands, result_filepaths

def parallel_blast_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles
    
        f: file containing one set of mapping instructions per line
        
        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log
         
        If f contained the two lines above, this function would 
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    lines = list(f)
    # handle catting of log files
    basic_process_run_results_f([lines[1]])
    # handle merging of otu maps
    fields = lines[0].strip().split()
    infiles_list = fields[:-1]
    out_filepath = fields[-1] 
    try:
        of = open(out_filepath,'w')
    except IOError:
        raise IOError,\
         "Poller can't open final output file: %s" % out_filepath  +\
         "\nLeaving individual jobs output.\n Do you have write access?"

    unique_otu_map = {}
    for fp in infiles_list:
        for line in open(fp):
            fields = line.strip().split()
            try:
                # current otu_id already exists, so append this
                # set of seq_ids
                unique_otu_map[fields[0]] += fields[1:]
            except KeyError:
                # current otu_id has not been seen yet, so 
                # create it with the current set of otus
                unique_otu_map[fields[0]] = fields[1:]
    
    for otu_id, seq_ids in unique_otu_map.items():
        of.write('\t'.join([otu_id] + seq_ids))
        of.write('\n')            
    of.close()
    
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback, 
    # you'll get an error right away rather than going into an infinite loop
    return True

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
