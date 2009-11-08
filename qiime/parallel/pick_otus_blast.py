#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# pick_otus_blast.py

""" Description
File created on 25 Aug 2009.

"""
from __future__ import division
from optparse import OptionParser
from os import popen, system, makedirs, mkdir
from os.path import split, splitext
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_rename_command, write_filepaths_to_file,\
    write_merge_map_file_pick_otus
from qiime.parallel.poller import basic_process_run_results_f
from qiime.util import qiime_config

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the Qiime Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

def get_job_commands(python_exe_fp,pick_otus_fp,fasta_fps,\
    output_dir,blast_db,job_prefix,working_dir,max_e_value,\
    command_prefix=None,command_suffix=None):
    """Generate pick_otus commands which should be submitted to cluster
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_otus.log', 
                     job_prefix + '.%d_otus.txt']
    
    # Initialize the command_prefix and command_suffix
    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
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
         '%s %s %s -i %s -b %s -m blast -o %s -e %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          pick_otus_fp,\
          fasta_fp,\
          blast_db,\
          working_dir,\
          max_e_value,\
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
    seconds_to_sleep,command_prefix=None,command_suffix=None):
    """Generate command to initiate a poller to monitior/process completed runs
    """
    
    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
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

usage_str = """usage: %prog [options] {-i INPUT_FP -o OUTPUT_DIR}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

Split the input file (-i) into five jobs (-O) to align against 
 pynast_test_template.fasta (-t), submit the jobs to the cluster (default)
 and write the output (-o) to /home/caporaso/out:

 python ~/Qiime/qiime/parallel/align_seqs_pynast.py -i 10_seq.fasta -O 5 -t /data/pynast_test_template.fasta -o /home/caporaso/out
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
          
    # define relevant align_seqs.py parameters
    parser.add_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp [REQUIRED]') 
    
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files '+\
           '[REQUIRED]')
           
    parser.add_option('-e','--max_e_value',\
          dest='max_e_value',help='Max E-value '+\
          '[default: %default]', default='1e-30')
           
    parser.add_option('-p','--pick_otus_fp',
            help='path to pick_otus.py [default: %default]',\
            default=qiime_config['pick_otus_fp'])
          
    parser.add_option('-r','--refseqs_fp',action='store',\
           type='string',help='full path to '+\
           'template alignment [default: %default]')
    
    parser.add_option('-b','--blast_db',action='store',\
           type='string',help='database to blast against '+\
           '[default: %default]')
          
    # Define parallel-script-specific parameters
    parser.add_option('-N','--align_seqs_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/align_seqs.py [default: %default]',\
           default=qiime_config['align_seqs_fp'])
        
    parser.add_option('-O','--jobs_to_start',type='int',\
            help='Number of jobs to start [default: %default]',default=24)
           
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
            
    parser.add_option('-U','--cluster_jobs_fp',
            help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default=qiime_config['cluster_jobs_fp'])

    parser.add_option('-W','--suppress_polling',action='store_true',
           help='suppress polling of jobs and merging of results '+\
           'upon completion [default: %default]',\
           default=False)
           
    parser.add_option('-X','--job_prefix',help='job prefix '+\
           '[default: ALIGN_ + 4 random chars]')

    parser.add_option('-Y','--python_exe_fp',
           help='full path to python executable [default: %default]',\
           default=qiime_config['python_exe_fp'])
        
    parser.add_option('-Z','--seconds_to_sleep',type='int',\
            help='Number of seconds to sleep between checks for run '+\
            ' completion when polling runs [default: %default]',default=60)
                             
    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp','output_dir']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
    
    if opts.blast_db == None and opts.refseqs_fp == None:
        parser.error('Either blast_db or refseqs_fp must be provided.')

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    pick_otus_fp = opts.pick_otus_fp
    refseqs_fp = opts.refseqs_fp
    cluster_jobs_fp = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    max_e_value = opts.max_e_value

    created_temp_paths = []

    if not opts.blast_db:        
        # Build the blast database from the reference_seqs_fp -- all procs
        # will then access one db rather than create one per proc
        blast_db, db_files_to_remove = \
             build_blast_db_from_fasta_path(refseqs_fp)
        created_temp_paths += db_files_to_remove
    else:
        blast_db = opts.blast_db
    
    # split the input filepath into directory and filename, base filename and
    # extension
    input_dir, input_fasta_fn = split(input_fasta_fp)
    input_file_basename, input_fasta_ext = splitext(input_fasta_fn)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with RDP
    job_prefix = opts.job_prefix or get_random_job_prefix('POTU')
    
    # A temporary output directory is created in output_dir named
    # job_prefix. Output files are then moved from the temporary 
    # directory to the output directory when they are complete, allowing
    # a poller to detect when runs complete by the presence of their
    # output files.
    working_dir = '%s/%s' % (output_dir,job_prefix)
    try:
        mkdir(working_dir)
        created_temp_paths.append(working_dir)
    except OSError:
        # working dir already exists
        pass
    
    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,\
      job_prefix,working_dir=output_dir)
    created_temp_paths += tmp_fasta_fps
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)
    
    # generate the list of commands to be pushed out to nodes and the list of
    # output files generated by each job
    commands, job_result_filepaths = \
     get_job_commands(python_exe_fp,pick_otus_fp,tmp_fasta_fps,\
     output_dir,blast_db,job_prefix,working_dir,max_e_value)
    created_temp_paths += job_result_filepaths

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
        process_run_results_f =\
         'qiime.parallel.pick_otus_blast.parallel_blast_process_run_results_f'
        write_merge_map_file_pick_otus(job_result_filepaths,output_dir,\
            merge_map_filepath,input_file_basename)
        created_temp_paths.append(merge_map_filepath)
        
        # Create the filepath listing the temporary files to be deleted,
        # but don't write it yet
        deletion_list_filepath = '%s/deletion_list.txt' % working_dir
        created_temp_paths.append(deletion_list_filepath)
        
        # Generate the command to run the poller, and the list of temp files
        # created by the poller
        poller_command, poller_result_filepaths =\
         get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
         merge_map_filepath,deletion_list_filepath,process_run_results_f,\
         seconds_to_sleep=seconds_to_sleep)
        created_temp_paths += poller_result_filepaths
        
        # append the poller command to the list of job commands
        commands.append(poller_command)
        
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
        submit_jobs(cluster_jobs_fp,jobs_fp,job_prefix)

    