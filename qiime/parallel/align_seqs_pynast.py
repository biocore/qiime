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
from qiime.util import qiime_config
from pynast.util import pairwise_alignment_methods

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the Qiime Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

def get_commands(python_exe_fp,align_seqs_fp,fasta_fps,template_aln_fp,\
    pairwise_alignment_method,output_dir,blast_db,blast_executable,\
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
    
    for fasta_fp in fasta_fps:
        command = \
         '%s %s %s %s -p %1.2f -e %d -b %s -m pynast -t %s -a %s -o %s -i %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          align_seqs_fp,\
          blast_option,\
          min_percent_id,\
          min_length,\
          blast_executable,\
          template_aln_fp,\
          pairwise_alignment_method,
          output_dir,
          fasta_fp,
          command_suffix)
          
        commands.append(command)
        
    return commands

usage_str = """usage: %prog [options] {-i INPUT_FP -t TEMPLATE_ALN_FP -o OUTPUT_DIR}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

Split the input file (-i) into five jobs (-j) to align against 
 pynast_test_template.fasta (-t), submit the jobs to the cluster (default)
 and write the output (-o) to /home/caporaso/out:

 parallel_pynast.py -i 10_seq.fasta -j 5 -t /data/pynast_test_template.fasta -o /home/caporaso/out
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
          
    parser.add_option('-t','--template_aln_fp',action='store',\
           type='string',help='full path to '+\
           'template alignment [REQUIRED]')
    
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files '+\
           '[REQUIRED]')
            
    pairwise_alignment_method_choices = pairwise_alignment_methods.keys()
    parser.add_option('-a','--pairwise_alignment_method',\
          type='choice',help='Method to use for pairwise alignments'+\
          ' (applicable with -m pynast) [default: %default]',\
          default='blast',choices=pairwise_alignment_method_choices)
    
    parser.add_option('-d','--blast_db',action='store',\
           type='string',help='database to blast against '+\
           '[default: formatted from template alignment]')
           
    parser.add_option('-b','--blastall_fp',
        default=qiime_config['blastall_fp'][0],
        help='Path to blastall [default: %default]')
        
    parser.add_option('-e','--min_length',action='store',\
          type='int',help='Minimum sequence '+\
          'length to include in alignment [default: %default]',\
          default=1000)
          
    parser.add_option('-p','--min_percent_id',action='store',\
          type='float',help='Minimum percent '+\
          'sequence identity to closest blast hit to include sequence in'+\
          ' alignment [default: %default]',default=75.0)
          
    # Define parallel-script-specific parameters
    parser.add_option('-n','--align_seqs_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/align_seqs.py [default: %default]',\
           default=qiime_config['align_seqs_fp'][0])
           
    parser.add_option('-f','--formatdb_fp',
        default=qiime_config['formatdb_fp'][0],
        help='Path to formatdb executable [default: %default]')
           
    parser.add_option('-x','--job_prefix',help='job prefix '+\
           '[default: ALIGN_ + 4 random chars]')
            
    parser.add_option('-u','--cluster_jobs_fp',
            help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default=qiime_config['cluster_jobs_fp'][0])

    parser.add_option('-S','--suppress_submit_jobs',action='store_true',\
            help='Only split input and write commands file - don\'t submit '+\
            'jobs [default: %default]',default=False)
        
    parser.add_option('-j','--jobs_to_start',type='int',\
            help='Number of jobs to start [default: %default]',default=24)

    parser.add_option('-y','--python_exe_fp',
           help='full path to python executable [default: %default]',\
           default=qiime_config['python_exe_fp'][0])
                             
    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp','template_aln_fp','output_dir']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    align_seqs_fp = opts.align_seqs_fp
    cluster_jobs_fp = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    template_aln_fp = opts.template_aln_fp
    pairwise_alignment_method = opts.pairwise_alignment_method
    formatdb_fp = opts.formatdb_fp
    blastall_fp = opts.blastall_fp
    min_length = opts.min_length
    min_percent_id = opts.min_percent_id
    

    blast_db = opts.blast_db
    ## If the user doesn't provide a blast database, one option is to
    ## to create it here so it doesn't get repeated on every proc. The 
    ## problem with creating it here is that it makes it harder to 
    ## clean-up because we don't know when all procs are done with it. 
    ## For now, will let each proc build the db itself.
    
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
    
    # generate the list of commands to be pushed out to nodes
    commands = \
     get_commands(python_exe_fp,align_seqs_fp,tmp_fasta_fps,template_aln_fp,\
     pairwise_alignment_method,output_dir,blast_db,blastall_fp,\
     min_length,min_percent_id)
     
    # write the commands to the 'jobs files'
    jobs_fp = write_jobs_file(commands,job_prefix=job_prefix)
    
    # submit the jobs file using cluster_jobs, if not suppressed by the
    # user
    if not opts.suppress_submit_jobs:
        submit_jobs(cluster_jobs_fp,jobs_fp,job_prefix)
    
    
    
    
    