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
__version__ = "1.2.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

def get_job_commands(python_exe_fp,pick_otus_fp,fasta_fps,
     output_dir,refseqs_fp,job_prefix,working_dir,similarity,
     enable_rev_strand_match,optimal_uclust,exact_uclust,
     max_accepts,max_rejects,stable_sort,save_uc_files,
     command_prefix='/bin/bash; ',
     command_suffix='; exit'):
    """Generate pick_otus commands which should be run
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_otus.log', 
                     job_prefix + '.%d_otus.txt',
                     job_prefix + '.%s_failures.txt']
    
    # Create lists to store the results
    commands = []
    result_filepaths = []
    
    if enable_rev_strand_match:
        enable_rev_strand_match_str = '-z'
    else:
        enable_rev_strand_match_str = ''
    if optimal_uclust:
        optimal_uclust_str = '-A'
    else:
        optimal_uclust_str = ''
    if exact_uclust:
        exact_uclust_str = '-E'
    else:
        exact_uclust_str = ''
    if stable_sort:
        stable_sort_str = '--uclust_stable_sort'
    else:
        stable_sort_str = ''
    if save_uc_files:
        save_uc_files = ''
        out_filenames += [job_prefix + '%d_clusters.uc']
    else:
        save_uc_files = '-d'
        
    
    # Iterate over the input files
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [fn % i for fn in out_filenames],working_dir,output_dir)
        result_filepaths += current_result_filepaths
            
        command = \
         '%s %s %s -i %s -r %s -m uclust_ref --suppress_new_clusters -o %s -s %s %s %s %s --max_accepts %s --max_rejects %s %s %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          pick_otus_fp,\
          fasta_fp,\
          refseqs_fp,\
          working_dir,\
          similarity,\
          enable_rev_strand_match_str,
          optimal_uclust_str,
          exact_uclust_str,
          max_accepts,
          max_rejects,
          stable_sort_str,
          save_uc_files,
          rename_command,
          command_suffix)

          
        commands.append(command)

    return commands, result_filepaths

def parallel_uclust_ref_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles
    
        f: file containing one set of mapping instructions per line
        
        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log
         f1_failures.txt f2_failures.txt f3_failures.txt f_failires.txt
         
        If f contained the two lines above, this function would 
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    lines = list(f)
    # handle catting of log files and failure files
    basic_process_run_results_f([lines[1],lines[2]])
    # # handle catting of failures files
    # basic_process_run_results_f([lines[2]])
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
