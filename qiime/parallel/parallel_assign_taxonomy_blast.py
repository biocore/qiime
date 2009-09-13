#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# parallel_assign_taxonomy_blast.py

""" Description
File created on 13 Sep 2009.

"""
from __future__ import division
from optparse import OptionParser
from os import popen, system
from os.path import split
from cogent.app.util import get_tmp_filename
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths

def get_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,reference_seqs_fp,\
    blastmat_path,fasta_fps,out_fps,log_fps,working_dir,\
    command_prefix=None,command_suffix=None):
    """Generate BlastTaxonAssiger classifier commands to be submitted to cluster
    """
    
    command_prefix = command_prefix or\
     '/bin/bash; cd %s; export BLASTMAT=%s;' \
       % (working_dir,blastmat_path)
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    
    for out_fp,log_fp,fasta_fp in zip(out_fps,log_fps,fasta_fps):
        command = '%s %s %s -o %s -l %s -m blast -r %s %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          assign_taxonomy_fp,\
          out_fp,
          log_fp,
          reference_seqs_fp,
          fasta_fp,
          id_to_taxonomy_fp,
          command_suffix)
          
        commands.append(command)
        
    return commands

usage_str = """usage: %prog [options] {-i INPUT_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

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

    parser.add_option('-u','--path_to_cluster_jobs',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default='/home/caporaso/bin/cluster_jobs.py')

    parser.add_option('-w','--working_dir',action='store',\
            type='string',help='directory to do work in' +\
            ' [default: %default]',\
            default='/home/caporaso/quicksand')


    parser.add_option('-r','--reference_seqs_fp',action='store',\
           type='string',help='full path to reference sequence filepath '+\
           '[default: %default]',\
           default='/quicksand/caporaso/silva_97/all_silva.fasta')
           
    parser.add_option('-a','--assign_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/assign_taxonomy.py [default: %default]',\
           default='/home/caporaso/Qiime/qiime/assign_taxonomy.py')
           
    parser.add_option('-t','--id_to_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'id_to_taxonomy mapping file [default: %default]',\
           default='/quicksand/caporaso/silva_97/silva_tax.txt')
           
    parser.add_option('-b','--blastmat_fp',action='store',\
           type='string',help='full path to '+\
           'blastmat file [default: %default]',\
           default='/home/caporaso/blast-2.2.16/data/')
                             
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
    path_to_cluster_jobs = opts.path_to_cluster_jobs
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    log_dir = opts.log_dir
    reference_seqs_fp = opts.reference_seqs_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    blastmat_fp = opts.blastmat_fp
    working_dir = opts.working_dir
    
    # split the input filepath into directory and filename
    input_dir, input_fasta_fn = split(input_fasta_fp)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with BTA (BlastTaxonAssigner)
    job_prefix = opts.job_prefix or get_random_job_prefix('BTA')
    
    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix,working_dir)
     
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
     get_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,reference_seqs_fp,\
     blastmat_fp,tmp_fasta_fps,out_fps,log_fps,working_dir)
     
    # write the commands to the 'jobs files'
    jobs_fp = write_jobs_file(commands,job_prefix=job_prefix)
    
    # submit the jobs file using cluster_jobs, if requested by the
    # user
    if opts.submit_jobs:
        submit_jobs(path_to_cluster_jobs,jobs_fp,job_prefix)
    
