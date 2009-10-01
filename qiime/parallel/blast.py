#!/usr/bin/env python
#parallel_blast.py: make and run parallel blast given file of seqs and db

__version__ = '0.01'

from os import mkdir, system
from os.path import exists
from optparse import OptionParser
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths
from qiime.util import qiime_config

def format_blast_database(db_path,formatdb_executable):
    system('%s -p F -i %s' % (formatdb_executable,db_path))

def get_commands(infile_paths,outfile_paths,db_path,blast_executable_path,\
    blastmat_path,e_value,word_size,num_hits,command_prefix=None,\
    command_suffix=None):
    
    command_prefix = command_prefix or\
     '/bin/bash; export BLASTMAT=%s;' % blastmat_path
    command_suffix = command_suffix or\
     '; exit'
     
    commands = []
    for infile_path, outfile_path in zip(infile_paths,outfile_paths):
        command = \
         "%s %s -p blastn -m 9 -e %s -W %s -b %s -i %s -d %s > %s %s" % \
         (command_prefix,\
          blast_executable_path,\
          e_value,\
          word_size,\
          num_hits, 
          infile_path,\
          db_path,\
          outfile_path,\
          command_suffix)
        commands.append(command)
    
    return commands
    
usage_str = """usage: %prog [options] {-i INPUT_FP -d  -d REFERENCE_SEQS_FP -o OUTPUT_DIR}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
 Split 10_seq.fasta (-i) into three fasta files (-j) and blast each against 
  blast database created from 1000_seq.fasta (-d):
 python Qiime/qiime/parallel/blast.py -i /home/caporaso/10_seq.fasta -d /home/caporaso/1000_seq.fasta -j 3 -o /home/caporaso/bla_out
"""
    
def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--infile_path',action='store',\
          type='string',dest='infile_path',
          help='Path of sequences to use as queries [REQUIRED]')

    parser.add_option('-d','--database_path',action='store',\
          type='string',dest='db_path',
            help='Path to database of sequences to search against' +\
            ' [REQUIRED]')
            
    parser.add_option('-o', '--output_dir', \
        help='name of output directory for blast jobs [REQUIRED]')
    
    parser.add_option('-e','--e_value',action='store',\
        type='float', default=1e-30, dest='e_value',
        help='E-value threshold for blasts [default: %default]')

    parser.add_option('-n','--num_hits',action='store',\
        type='int', default=1, dest='num_hits',
        help='number of hits per query for blast results [default: %default]')

    parser.add_option('-W','--word_size',action='store',\
        type='int', default=30, dest='word_size',
        help='word size for blast searches [default: %default]')

    parser.add_option('-D', '--suppress_format_blastdb', action='store_true',\
        default=False,help='supress format of blastdb [default: %default]')
    
    parser.add_option('-a','--blastmat_dir',action='store',\
           type='string',help='full path to directory containing '+\
           'blastmat file [default: %default]',\
           default=qiime_config['blastmat_dir'])   
           
    parser.add_option('-f','--formatdb_fp',
        default=qiime_config['formatdb_fp'],
        help='Path to formatdb executable [default: %default]')
           
    parser.add_option('-b','--blastall_fp',
        default=qiime_config['blastall_fp'],
        help='Path to blastall [default: %default]')

    # Define parallel-script-specific parameters
    parser.add_option('-x', '--job_prefix', action='store', \
        type='string', dest='job_prefix', default=None,
        help='job prefix, max 10 chars [default: %default]')
    
    parser.add_option('-u','--cluster_jobs_fp',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default=qiime_config['cluster_jobs_fp'])

    parser.add_option('-S','--suppress_submit_jobs',action='store_true',\
            help='Only split input and write commands file - don\'t submit '+\
            'jobs [default: %default]',default=False)
    
    parser.add_option('-j','--jobs_to_start',action='store',type='int',\
            help='Number of jobs to start [default: %default]',default=24)

    opts,args = parser.parse_args()
    
    required_options = ['infile_path','db_path','output_dir']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
            
    return opts, args


if __name__ == '__main__':
    from sys import exit
    opts, args = parse_command_line_parameters()
    
    # create local copies of command-line options
    input_fasta_fp = opts.infile_path 
    db_path = opts.db_path
    num_jobs_to_start = opts.jobs_to_start
    blastall_fp = opts.blastall_fp
    blastmat_dir = opts.blastmat_dir
    cluster_jobs_fp = opts.cluster_jobs_fp
    e_value = opts.e_value
    word_size = opts.word_size
    num_hits = opts.num_hits
    output_dir = opts.output_dir
    formatdb_fp = opts.formatdb_fp
    suppress_format_blastdb = opts.suppress_format_blastdb
    
    job_prefix = opts.job_prefix or get_random_job_prefix('BLAST')
    
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,num_jobs_to_start)
    
    if not suppress_format_blastdb:
        blast_db, db_files_to_remove = \
         build_blast_db_from_fasta_path(db_path)
        
    # split the fasta files and get the list of resulting files        
    tmp_fasta_fps =\
     split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix,output_dir)
      
    # build the list of output filepaths from the set of input files
    # by appending '.blast_out.txt' to the end of each
    out_fps = build_filepaths_from_filepaths(tmp_fasta_fps,\
     directory=output_dir,suffix='.blast_out.txt')

    # generate the list of commands to be pushed out to nodes    
    commands = get_commands(tmp_fasta_fps,out_fps,db_path,\
     blastall_fp,blastmat_dir,e_value,word_size,num_hits,\
     command_prefix=None,command_suffix=None)
    
    # write the commands to the 'jobs files'
    jobs_fp = write_jobs_file(commands,job_prefix=job_prefix)
    
    # submit the jobs file using cluster_jobs, if not suppressed by the
    # user
    if not opts.suppress_submit_jobs:
        submit_jobs(cluster_jobs_fp,jobs_fp,job_prefix)
