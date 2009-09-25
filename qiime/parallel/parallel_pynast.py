#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# parallel_pynast.py

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
from pynast.util import build_temp_blast_db_from_alignment_fp

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the Qiime Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

def get_commands(python_exe_fp,align_seqs_fp,fasta_fps,template_aln_fp,\
    pairwise_alignment_method,out_fps,log_fps,blast_db,blast_executable,\
     min_length,min_percent_id,command_prefix=None,command_suffix=None):
    """Generate PyNAST commands which should be submitted to cluster
    """
    
    command_prefix = command_prefix or '/bin/bash; '
    command_suffix = command_suffix or '; exit'
    
    commands = []
    
    # If there is a value for blast_db, pass it. If not, it
    # will be created on-the-fly. Note that on-the-fly blast dbs
    # are created with a string of random chars in the name, so this is safe.
    # They shouldn't overwrite one another, and will be cleaned up.
    if blast_db:
        blast_option = '-d %s' % blast_db
    else:
        blast_option = ''
    
    for out_fp,log_fp,fasta_fp in zip(out_fps,log_fps,fasta_fps):
        command = \
         '%s %s %s %s -p %1.2f -e %d -b %s -m pynast -t %s -a %s -o %s -l %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          align_seqs_fp,\
          blast_option,\
          min_percent_id,\
          min_length,\
          blast_executable,\
          template_aln_fp,\
          pairwise_alignment_method,
          out_fp,
          log_fp,
          fasta_fp,
          command_suffix)
          
        commands.append(command)
        
    return commands

usage_str = """usage: %prog [options] {-i INPUT_FP -t TEMPLATE_ALN_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

Split the input file (-i) into five jobs (-j) to align against 
 pynast_test_template.fasta (-t) and submit the jobs (-s)
 to the cluster:

 parallel_pynast.py -s -i 10_seq.fasta -j 5 -t /data/pynast_test_template.fasta
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)
          
    parser.add_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp [REQUIRED]') 
          
    parser.add_option('-t','--template_aln_fp',action='store',\
           type='string',help='full path to '+\
           'template alignment [REQUIRED]') 

    parser.add_option('-s','--submit_jobs',action='store_true',\
            help='Submit jobs in addition to splitting input and '+\
            'writing commands file  [default: %default]',default=False)
        
    parser.add_option('-j','--jobs_to_start',action='store',type='int',\
            help='Number of jobs to start [default: %default]',default=24)

    parser.add_option('-y','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default='/home/caporaso/bin/python')
    
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files '+\
           '[default: %default]',default='.')
    
    parser.add_option('-d','--blast_db',action='store',\
           type='string',help='database to blast against '+\
           '[default: formatted from template alignment]')
           
    parser.add_option('-b','--blast_executable',action='store',\
        type='string',
        default='/home/caporaso/bin/blastall',
        help='Path to blast executable [default: %default]')
        
    parser.add_option('-f','--formatdb_executable',action='store',\
        type='string',
        default='/home/caporaso/bin/formatdb',
        help='Path to formatdb executable [default: %default]')
           
    parser.add_option('-x','--job_prefix',action='store',\
           type='string',help='job prefix '+\
           '[default: ALIGN_ + 4 random chars]')
    
    parser.add_option('-l','--log_dir',action='store',\
           type='string',help='path to store log files '+\
           '[default: %default]',default='.')
           
    parser.add_option('-n','--align_seqs_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/assign_taxonomy.py [default: %default]',\
           default='/home/caporaso/Qiime/qiime/align_seqs.py')
          
    parser.add_option('-a','--pairwise_alignment_method',action='store',\
          type='string',help='Method to use for pairwise alignments'+\
          ' (applicable with -m pynast) [default: %default]',\
          default='blast')
            
    parser.add_option('-u','--path_to_cluster_jobs',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default='/home/caporaso/bin/cluster_jobs.py')
            
    parser.add_option('-e','--min_length',action='store',\
          type='int',help='Minimum sequence '+\
          'length to include in alignment [default: %default]',\
          default=1000)
          
    parser.add_option('-p','--min_percent_id',action='store',\
          type='float',help='Minimum percent '+\
          'sequence identity to closest blast hit to include sequence in'+\
          ' alignment [default: %default]',default=75.0)
                             
    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp','template_aln_fp']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    align_seqs_fp = opts.align_seqs_fp
    path_to_cluster_jobs = opts.path_to_cluster_jobs
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    log_dir = opts.log_dir
    template_aln_fp = opts.template_aln_fp
    pairwise_alignment_method = opts.pairwise_alignment_method
    formatdb_executable = opts.formatdb_executable
    blast_executable = opts.blast_executable
    min_length = opts.min_length
    min_percent_id = opts.min_percent_id
    

    blast_db = opts.blast_db
    ## If the user doesn't provide a blast database, one option is to
    ## to create it here so it doesn't get repeated on every proc. The 
    ## problem with creating it here is that it makes it harder to 
    ## clean-up. For now, will let each proc build the db itself.
    # if not blast_db:
    #     blast_db, files_to_remove = \
    #       build_temp_blast_db_from_alignment_fp(\
    #       template_aln_fp,formatdb_executable)
    #     blast_db = str(blast_db)
    
    # split the input filepath into directory and filename
    input_dir, input_fasta_fn = split(input_fasta_fp)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with RDP
    job_prefix = opts.job_prefix or get_random_job_prefix('ALIGN')
    
    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix)
     
    # build the list of output filepaths from the set of input files
    # by appending '.tax_assignments.txt' to the end of each
    out_fps = build_filepaths_from_filepaths(tmp_fasta_fps,\
     directory=output_dir,suffix='.aligned.fasta')
     
    # build the list of log filepaths from the set of input files
    # by appending '.tax_assignments.log' to the end of each
    log_fps = build_filepaths_from_filepaths(tmp_fasta_fps,\
     directory=log_dir,suffix='.aligned.log')
    
    # generate the list of commands to be pushed out to nodes
    commands = \
     get_commands(python_exe_fp,align_seqs_fp,tmp_fasta_fps,template_aln_fp,\
     pairwise_alignment_method,out_fps,log_fps,blast_db,blast_executable,\
     min_length,min_percent_id)
     
    # write the commands to the 'jobs files'
    jobs_fp = write_jobs_file(commands,job_prefix=job_prefix)
    
    # submit the jobs file using cluster_jobs, if requested by the
    # user
    if opts.submit_jobs:
        submit_jobs(path_to_cluster_jobs,jobs_fp,job_prefix)
    
    
    
    
    