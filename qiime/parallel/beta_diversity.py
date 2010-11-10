#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# beta_diversity.py

""" Description
File created on 7 Jan 2010.

"""
from __future__ import division
from os.path import split, join
import os
from qiime.parallel.util import get_rename_command, merge_to_n_commands
from qiime.parse import parse_otu_table
from qiime.format import format_distance_matrix

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso","Justin Kuczynski"] 
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

def get_job_commands_multiple_otu_tables(
    python_exe_fp,beta_diversity_fp,tree_fp,job_prefix,metrics,input_fps,
    output_dir,working_dir,command_prefix=None,command_suffix=None,
    full_tree=False):
    """Generate beta diversity to split multiple OTU tables to multiple jobs
    """

    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    if full_tree:
        full_tree_str = '-f'
    else:
        full_tree_str = ''
    
    commands = []
    result_filepaths = []
    
    for input_fp in input_fps:
        input_path, input_fn = split(input_fp)
        output_fns = ['%s_%s' % (metric, input_fn) \
         for metric in metrics.split(',')]
        rename_command, current_result_filepaths = get_rename_command(\
         output_fns,working_dir,output_dir)
        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -i %s -o %s -t %s -m %s %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          beta_diversity_fp,\
          input_fp,
          working_dir + '/',
          tree_fp,
          metrics,
          full_tree_str,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths
    
def get_job_commands_single_otu_table(
    python_exe_fp,beta_diversity_fp,tree_fp,job_prefix,metrics,input_fp,
    output_dir,working_dir,jobs_to_start,command_prefix=None,
    command_suffix=None):
    """Generate beta diversity to split single OTU table to multiple jobs
    
    always passes -f to beta_diversity.py
    """

    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    commands = []
    result_filepaths = []
    
    sids = parse_otu_table(open(input_fp,'U'))[0]
    
    sample_id_groups = merge_to_n_commands(sids,jobs_to_start,',','','')
    for i, sample_id_group in enumerate(sample_id_groups):
        working_dir_i = os.path.join(working_dir, str(i))
        output_dir_i = os.path.join(output_dir, str(i))
        input_dir, input_fn = split(input_fp)
        sample_id_desc = sample_id_group.replace(',','_')
        output_fns = ['%s_%s' % (metric, input_fn) \
         for metric in metrics.split(',')]
        rename_command, current_result_filepaths = get_rename_command(\
         output_fns,working_dir_i,output_dir_i)

        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -i %s -o %s -t %s -m %s -f -r %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          beta_diversity_fp,\
          input_fp,
          working_dir_i + '/',
          tree_fp,
          metrics,
          sample_id_group,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths
    
def create_merge_map_file_single_otu_table(input_fp,output_dir,
    metrics,merge_map_filepath,expected_files_filepath):
    
    merge_map_f = open(merge_map_filepath,'w')
    
    input_dir, input_fn = split(input_fp)
    
    expected_output_files = [fp.strip() for fp in 
     open(expected_files_filepath,'U')]
    
    for metric in metrics.split(','):
        fps_to_merge = [fp for fp in expected_output_files if '/%s_' % metric in fp]
        output_fp = join(output_dir,'%s_%s' % (metric,input_fn))
        merge_map_f.write('%s\t%s\n' % ('\t'.join(fps_to_merge),output_fp))
    
    merge_map_f.close()

def assemble_distance_matrix(dm_components):
    """ assemble distance matrix components into a complete dm string
    
    """
    data = {}
    # iterate over compenents
    for c in dm_components:
        # create a blank list to store the column ids
        col_ids = []
        # iterate over lines
        for line in c:
            # split on tabs remove leading and trailing whitespace
            fields = line.strip().split()
            if fields:
                # if no column ids seen yet, these are them
                if not col_ids:
                    col_ids = fields
                # otherwise this is a data row so add it to data
                else:
                    sid = fields[0]
                    data[sid] = dict(zip(col_ids,fields[1:]))

    # grab the col/row ids as a list so it's ordered
    labels = data.keys()
    # create an empty list to build the dm
    dm = []
    # construct the dm one row at a time
    for l1 in labels:
        dm.append([data[l1][l2] for l2 in labels])
    # create the dm string and return it
    dm = format_distance_matrix(labels,dm)
    return dm
        
def parallel_beta_diversity_process_run_results_f(f):
    """ Handles re-assembling of an OTU table from component vectors
    """
    # iterate over component, output fp lines
    for line in f:
        fields = line.strip().split('\t')
        dm_components = fields[:-1]
        output_fp = fields[-1]
        # assemble the current dm
        dm = assemble_distance_matrix(map(open,dm_components))
        # and write it to file
        output_f = open(output_fp,'w')
        output_f.write(dm)
        output_f.close()
        
    return True

def get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
    merge_map_filepath,deletion_list_filepath,seconds_to_sleep,\
    process_run_results_f,
    command_prefix='/bin/bash; ',command_suffix='; exit'):
    """Generate command to initiate a poller to monitior/process completed runs
    """
    
    if process_run_results_f != None:
        process_run_results_f_str = '-p %s' % process_run_results_f
    else:
        process_run_results_f_str = ''
    
    result = '%s %s %s -f %s -m %s -d %s -t %d %s %s' % \
     (command_prefix,
      python_exe_fp,
      poller_fp,
      expected_files_filepath,
      merge_map_filepath,
      deletion_list_filepath,
      seconds_to_sleep,
      process_run_results_f_str,
      command_suffix)
      
    return result, []
