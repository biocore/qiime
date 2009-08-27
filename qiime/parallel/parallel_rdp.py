#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# parallel_rdp.py

""" Description
File created on 25 Aug 2009.

"""
from __future__ import division
from optparse import OptionParser
from os import popen, system
from os.path import split
from cogent.app.util import get_tmp_filename
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths

def get_commands(python_exe_fp,assign_taxonomy_fp,confidence,fasta_fps,\
    rdp_jar_fp,out_fps,log_fps,command_prefix=None,command_suffix=None):
    """Generate RDP classifier commands which should be submitted to cluster
    """
    
    command_prefix = command_prefix or\
     '/bin/bash; export RDP_JAR_PATH=%s; ' % rdp_jar_fp
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    
    for out_fp,log_fp,fasta_fp in zip(out_fps,log_fps,fasta_fps):
        command = '%s %s %s -c %1.2f -m rdp -o %s -l %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          assign_taxonomy_fp,\
          confidence,
          out_fp,
          log_fp,
          fasta_fp,
          command_suffix)
          
        commands.append(command)
        
    return commands

usage_str = """usage: %prog [options] {-i INPUT_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

Split the input file (-i) into five jobs (-j) and submit them (-s)
 to the cluster:

 parallel_rdp.py -j 5 -i /home/caporaso/10_seq.fasta -s
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)
          
    parser.add_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp [REQUIRED]') 

    parser.add_option('-s','--submit_jobs',action='store_true',\
            help='Submit jobs in addition to splitting input and '+\
            'writing commands file  [default: %default]',default=False)
        
    parser.add_option('-j','--jobs_to_start',action='store',type='int',\
            help='Number of jobs to start [default: %default]',default=24)

    parser.add_option('-p','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default='/home/caporaso/bin/python')
    
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files '+\
           '[default: %default]',default='.')
           
    parser.add_option('-x','--job_prefix',action='store',\
           type='string',help='job prefix '+\
           '[default: RDP_ + 5 random chars]')
    
    parser.add_option('-l','--log_dir',action='store',\
           type='string',help='path to store log files '+\
           '[default: %default]',default='.')

    parser.add_option('-r','--rdp_classifier_jar_path',action='store',\
           type='string',help='full path to rdp classifier jar file '+\
           '[default: %default]',\
           default='/home/caporaso/rdp_classifier/rdp_classifier-2.0.jar')
           
    parser.add_option('-a','--assign_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/assign_taxonomy.py [default: %default]',\
           default='/home/caporaso/Qiime/qiime/assign_taxonomy.py')
          
    parser.add_option('-c','--confidence',action='store',\
          type='float',help='Minimum confidence to'+\
          ' record an assignment [default: %default]',default=0.80)
            
    parser.add_option('-u','--path_to_cluster_jobs',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default='/home/caporaso/bin/cluster_jobs.py')
                             
    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    assign_taxonomy_fp = opts.assign_taxonomy_fp
    confidence = opts.confidence
    rdp_classifier_jar_path = opts.rdp_classifier_jar_path
    path_to_cluster_jobs = opts.path_to_cluster_jobs
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    log_dir = opts.log_dir
    
    # split the input filepath into directory and filename
    input_dir, input_fasta_fn = split(input_fasta_fp)
    
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with RDP
    job_prefix = opts.job_prefix or get_random_job_prefix('RDP')
    
    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix)
     
    # build the list of output filepaths from the set of input files
    # by appending '.tax_assignments.txt' to the end of each
    out_fps = build_filepaths_from_filepaths(tmp_fasta_fps,\
     directory=output_dir,suffix='.tax_assignments.txt')
     
    # build the list of log filepaths from the set of input files
    # by appending '.tax_assignments.log' to the end of each
    log_fps = build_filepaths_from_filepaths(tmp_fasta_fps,\
     directory=log_dir,suffix='.tax_assignments.log')
    
    # generate the list of commands to be pushed out to nodes
    commands = \
     get_commands(python_exe_fp,assign_taxonomy_fp,confidence,tmp_fasta_fps,\
     rdp_classifier_jar_path,out_fps,log_fps)
     
    # write the commands to the 'jobs files'
    jobs_fp = write_jobs_file(commands,job_prefix=job_prefix)
    
    # submit the jobs file using cluster_jobs, if requested by the
    # user
    if opts.submit_jobs:
        submit_jobs(path_to_cluster_jobs,jobs_fp,job_prefix)
    
    
    
    
    