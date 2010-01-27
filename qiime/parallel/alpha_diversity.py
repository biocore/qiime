#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# alpha_diversity.py

""" Description
File created on 13 Sep 2009.

"""
from __future__ import division
from optparse import OptionParser
from glob import glob
from os import popen, system, mkdir, makedirs
from os.path import split, splitext
from subprocess import check_call, CalledProcessError
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_poller_command, get_rename_command, write_filepaths_to_file,\
    write_merge_map_file_assign_taxonomy
from qiime.alpha_diversity import list_known_metrics
from qiime.util import qiime_config

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso","Justin Kuczynski"] 
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

def get_job_commands(python_exe_fp,alpha_diversity_fp,tree_fp,job_prefix,\
    metrics,input_fps,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    """Generate alpha diversity commands to be submitted to cluster
    """

    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    commands = []
    result_filepaths = []
    
    for input_fp in input_fps:
        input_path, input_fn = split(input_fp)
        output_fn = 'alpha_%s' % input_fn
        rename_command, current_result_filepaths = get_rename_command(\
         [output_fn],working_dir,output_dir)
        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -i %s -o %s -t %s -m %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          alpha_diversity_fp,\
          input_fp,
          working_dir + '/' + output_fn,
          tree_fp,
          metrics,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths

usage_str = """usage: %prog [options] [{-i INPUT_FILE -o OUTPUT_PATH -t TREE_PATH -M METRICS} | -s]

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

List available metrics:
python alpha_diversity.py -s

Apply the observed_species, chao1, PD_whole_tree metrics (-m) to all otu tables 
 in ./rare/ (-i) and write the resulting output files to ./out/ (-o, will be 
 created if it doesn't exist). Use the tree file rep_set.tre (-t) when necessary. 

 python alpha_diversity.py -i ./rare/ -o ./out -m observed_species,chao1,PD_whole_tree -t ./rep_set.tre
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)
    
    # define relevant alpha_diversity.py parameters
    parser.add_option('-i', '--input_path',
        help='input path, must be directory [REQUIRED]')
        
    parser.add_option('-o', '--output_path',
        help='output path, must be directory [REQUIRED]')

    parser.add_option('-m', '--metrics',
        help='metrics to use, comma delimited [REQUIRED]')

    parser.add_option('-s', '--show_metrics', action='store_true', 
        dest="show_metrics",
        help='show available alpha diversity metrics and quit')
    
    parser.add_option('-t', '--tree_path',
        help='path to newick tree file, required for phylogenetic metrics'+\
        ' [default: %default]')   
           
    # Define parallel-script-specific parameters
    parser.add_option('-N','--alpha_diversity_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/alpha_diversity.py [default: %default]',\
           default=qiime_config['alpha_diversity_fp'])
        
    # Currently don't have control over number of jobs -- this is determined
    # by the number of rarefied data sets.
    # parser.add_option('-O','--jobs_to_start',action='store',type='int',\
    #         help='Number of jobs to start [default: %default]',default=24)
           
    parser.add_option('-P','--poller_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/parallel/poller.py [default: %default]',\
           default=qiime_config['poller_fp'])
           
    parser.add_option('-R','--retain_temp_files',action='store_true',\
           help='retain temporary files after runs complete '+\
           '(useful for debugging) [default: %default]',\
           default=False)
           
    parser.add_option('-S','--suppress_submit_jobs',action='store_true',\
            help='Only split input and write commands file - don\'t submit '+\
            'jobs [default: %default]',default=False)

    parser.add_option('-T','--poll_directly',action='store_true',\
            help='Poll directly for job completion rather than running '+\
            'poller as a separate job. If -T is specified this script will '+\
            'not return until all jobs have completed. [default: %default]',\
            default=False)

    parser.add_option('-U','--cluster_jobs_fp',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default=qiime_config['cluster_jobs_fp'])

    parser.add_option('-W','--suppress_polling',action='store_true',
           help='suppress polling of jobs and merging of results '+\
           'upon completion [default: %default]',\
           default=False)
           
    parser.add_option('-X','--job_prefix',action='store',\
           type='string',help='job prefix '+\
           '[default: ALDIV_ + 3 random chars]')
           
    parser.add_option('-Y','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default=qiime_config['python_exe_fp'])
        
    parser.add_option('-Z','--seconds_to_sleep',type='int',\
            help='Number of seconds to sleep between checks for run '+\
            ' completion when polling runs [default: %default]',default=60)
                             
    opts,args = parser.parse_args()
    
    if opts.show_metrics:
        print("Known metrics are: %s\n" \
                % (', '.join(list_known_metrics()),))
        exit(0)
    
    required_options = ['input_path','output_path','metrics']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of command-line options
    input_dir = opts.input_path
    output_dir = opts.output_path
    metrics = opts.metrics
    tree_fp = opts.tree_path
    
    alpha_diversity_fp = opts.alpha_diversity_fp
    python_exe_fp = opts.python_exe_fp
    path_to_cluster_jobs = opts.cluster_jobs_fp
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    poll_directly = opts.poll_directly

    created_temp_paths = []
    input_fps = glob('%s/*' % input_dir)
    # split the input filepath into directory and filename, base filename and
    # extension
    # input_dir, input_fn = split(input_path)
    # input_file_basename, input_file_ext = splitext(input_fn)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with ALDIV (ALphaDIVersity)
    job_prefix = opts.job_prefix or get_random_job_prefix('ALDIV')
    
    # A temporary output directory is created in output_dir named
    # job_prefix. Output files are then moved from the temporary 
    # directory to the output directory when they are complete, allowing
    # a poller to detect when runs complete by the presence of their
    # output files.
    working_dir = '%s/%s' % (output_dir,job_prefix)
    try:
        makedirs(working_dir)
        created_temp_paths.append(working_dir)
    except OSError:
        # working dir already exists
        pass
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)
    
    # Get the list of commands to be run and the expected result files
    commands, job_result_filepaths  = \
     get_job_commands(python_exe_fp,alpha_diversity_fp,tree_fp,job_prefix,\
     metrics,input_fps,output_dir,working_dir)
    
    # Set up poller apparatus if the user does not suppress polling
    if not suppress_polling:
        # Write the list of files which must exist for the jobs to be 
        # considered complete
        expected_files_filepath = '%s/expected_out_files.txt' % working_dir
        write_filepaths_to_file(job_result_filepaths,expected_files_filepath)
        created_temp_paths.append(expected_files_filepath)
        
        # Write the mapping file which described how the output files from
        # each job should be merged into the final output files
        merge_map_filepath = '%s/merge_map.txt' % working_dir
        open(merge_map_filepath,'w').close()
        created_temp_paths.append(merge_map_filepath)
        
        # Create the filepath listing the temporary files to be deleted,
        # but don't write it yet
        deletion_list_filepath = '%s/deletion_list.txt' % working_dir
        created_temp_paths.append(deletion_list_filepath)
        
        # Generate the command to run the poller, and the list of temp files
        # created by the poller
        if not poll_directly:
            poller_command, poller_result_filepaths =\
             get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
             merge_map_filepath,deletion_list_filepath,\
             seconds_to_sleep=seconds_to_sleep)
            # append the poller command to the list of job commands
            commands.append(poller_command)
        else:
            poller_command, poller_result_filepaths =\
             get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
             merge_map_filepath,deletion_list_filepath,\
             seconds_to_sleep=seconds_to_sleep,\
             command_prefix='',command_suffix='')
            
        created_temp_paths += poller_result_filepaths
        
        if not retain_temp_files:
            # If the user wants temp files deleted, now write the list of 
            # temp files to be deleted
            write_filepaths_to_file(created_temp_paths,deletion_list_filepath)
        else:
            # Otherwise just write an empty file
            write_filepaths_to_file([],deletion_list_filepath)
    
    # write the commands to the 'jobs files'
    write_jobs_file(commands,job_prefix=job_prefix,jobs_fp=jobs_fp)
    
    # submit the jobs file using cluster_jobs, if not suppressed by the
    # user
    if not opts.suppress_submit_jobs:
        submit_jobs(path_to_cluster_jobs,jobs_fp,job_prefix)

    if poll_directly:
        try:
            check_call(poller_command.split())
        except CalledProcessError, e:
            print '**Error occuring when calling the poller directly. '+\
            'Jobs may have been submitted, but are not being polled.'
            print str(e)
            exit(-1)
