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
from qiime.util import qiime_config

def get_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,reference_seqs_fp,\
    blastmat_path,fasta_fps,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    """Generate BlastTaxonAssiger classifier commands to be submitted to cluster
    """
    
    command_prefix = command_prefix or\
     '/bin/bash; cd %s; export BLASTMAT=%s;' \
       % (working_dir,blastmat_path)
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    
    for fasta_fp in fasta_fps:
        command = '%s %s %s -o %s -m blast -r %s -i %s -t %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          assign_taxonomy_fp,\
          output_dir,
          reference_seqs_fp,
          fasta_fp,
          id_to_taxonomy_fp,
          command_suffix)
          
        commands.append(command)
        
    return commands

usage_str = """usage: %prog [options] {-i INPUT_FP -o OUTPUT_DIR -t ID_TO_TAXONOMY_FP -r REFERENCE_SEQS_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
 python Qiime/qiime/parallel/parallel_assign_taxonomy_blast.py -j 5 -i /home/caporaso/at_inseqs.fasta -t /home/caporaso/at_id_to_taxonomy.txt -r /home/caporaso/at_refseqs.fasta -o /home/caporaso/out/


"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)
    
    # define relevant align_seqs.py parameters
    parser.add_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp [REQUIRED]') 
           
    parser.add_option('-t','--id_to_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'id_to_taxonomy mapping file [REQUIRED]')

    parser.add_option('-r','--reference_seqs_fp',action='store',\
           type='string',help='full path to reference sequence filepath '+\
           '[REQUIRED]')
    
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files '+\
           '[REQUIRED]')
           
    # Define parallel-script-specific parameters
    parser.add_option('-S','--suppress_submit_jobs',action='store_true',\
            help='Only split input and write commands file - don\'t submit '+\
            'jobs [default: %default]',default=False)
        
    parser.add_option('-j','--jobs_to_start',action='store',type='int',\
            help='Number of jobs to start [default: %default]',default=24)

    parser.add_option('-p','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default=qiime_config['python_exe_fp'])
           
    parser.add_option('-x','--job_prefix',action='store',\
           type='string',help='job prefix '+\
           '[default: BTA_ + 5 random chars]')

    parser.add_option('-u','--cluster_jobs_fp',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default=qiime_config['cluster_jobs_fp'])

    parser.add_option('-w','--working_dir',action='store',\
            type='string',help='directory to do work in' +\
            ' [default: %default]',\
            default=qiime_config['working_dir'])
           
    parser.add_option('-a','--assign_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/assign_taxonomy.py [default: %default]',\
           default=qiime_config['assign_taxonomy_fp'])
           
    parser.add_option('-b','--blastmat_dir',action='store',\
           type='string',help='full path to directory containing '+\
           'blastmat file [default: %default]',\
           default=qiime_config['blastmat_dir'])
                             
    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp','id_to_taxonomy_fp',\
     'reference_seqs_fp','output_dir']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    assign_taxonomy_fp = opts.assign_taxonomy_fp
    path_to_cluster_jobs = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    reference_seqs_fp = opts.reference_seqs_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    blastmat_fp = opts.blastmat_dir
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
    
    # generate the list of commands to be pushed out to nodes
    commands = \
     get_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,reference_seqs_fp,\
     blastmat_fp,tmp_fasta_fps,output_dir,working_dir)
     
    # write the commands to the 'jobs files'
    jobs_fp = write_jobs_file(commands,job_prefix=job_prefix)
    
    # submit the jobs file using cluster_jobs, if not suppressed by the
    # user
    if not opts.suppress_submit_jobs:
        submit_jobs(path_to_cluster_jobs,jobs_fp,job_prefix)
    
