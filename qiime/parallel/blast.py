#!/usr/bin/env python
#parallel_blast.py: make and run parallel blast given file of seqs and db

__author__ = "Rob Knight, Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from os import makedirs, system
from os.path import exists, split, splitext
from subprocess import check_call, CalledProcessError
from optparse import OptionParser
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_poller_command, get_rename_command, write_filepaths_to_file,\
    write_merge_map_file_blast
from qiime.util import qiime_config

def get_commands(infile_paths,db_path,blast_executable_path,\
    blastmat_path,e_value,word_size,num_hits,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    
    command_prefix = command_prefix or\
     '/bin/bash; export BLASTMAT=%s;' % blastmat_path
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    result_filepaths = []
    
    for i, infile_path in enumerate(infile_paths):
        
        infile_basename = splitext(split(infile_path)[1])[0]
        working_outfile_path = '%s/%s_blast_out.txt' %\
          (working_dir,infile_basename)
        outfile_path = '%s/%s_blast_out.txt' % (output_dir,infile_basename)
        
        rename_command = '; mv %s %s' % (working_outfile_path, outfile_path)
        
        result_filepaths.append(outfile_path)
        
        command = \
         "%s %s -p blastn -m 9 -e %s -W %s -b %s -i %s -d %s > %s %s %s" % \
         (command_prefix,\
          blast_executable_path,\
          e_value,\
          word_size,\
          num_hits, 
          infile_path,\
          db_path,\
          working_outfile_path,\
          rename_command,\
          command_suffix)
        commands.append(command)
    
    return commands, result_filepaths
    
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

    parser.add_option('-w','--word_size',action='store',\
        type='int', default=30, dest='word_size',
        help='word size for blast searches [default: %default]')

    parser.add_option('-D', '--suppress_format_blastdb', action='store_true',\
        default=False,help='supress format of blastdb [default: %default]')
    
    parser.add_option('-a','--blastmat_dir',action='store',\
           type='string',help='full path to directory containing '+\
           'blastmat file [default: %default]',\
           default=qiime_config['blastmat_dir'])   
           
    parser.add_option('-b','--blastall_fp',
        default=qiime_config['blastall_fp'],
        help='Path to blastall [default: %default]')

    # Define parallel-script-specific parameters
    
    parser.add_option('-O','--jobs_to_start',action='store',type='int',\
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
           
    parser.add_option('-X', '--job_prefix', action='store', \
        type='string', dest='job_prefix', default=None,
        help='job prefix, max 10 chars [default: %default]')
            
    parser.add_option('-Y','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default=qiime_config['python_exe_fp'])
        
    parser.add_option('-Z','--seconds_to_sleep',type='int',\
            help='Number of seconds to sleep between checks for run '+\
            ' completion when polling runs [default: %default]',default=60)

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
    python_exe_fp = opts.python_exe_fp
    num_jobs_to_start = opts.jobs_to_start
    blastall_fp = opts.blastall_fp
    blastmat_dir = opts.blastmat_dir
    cluster_jobs_fp = opts.cluster_jobs_fp
    e_value = opts.e_value
    word_size = opts.word_size
    num_hits = opts.num_hits
    output_dir = opts.output_dir
    suppress_format_blastdb = opts.suppress_format_blastdb
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    poll_directly = opts.poll_directly

    created_temp_paths = []
    
    # split the input filepath into directory and filename, base filename and
    # extension
    input_dir, input_fasta_fn = split(input_fasta_fp)
    input_file_basename, input_fasta_ext = splitext(input_fasta_fn)

    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with BLAST
    job_prefix = opts.job_prefix or get_random_job_prefix('BLAST')
    
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
        
    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,num_jobs_to_start)
    
    # Build the blast database if necessary
    if not suppress_format_blastdb:
        blast_db, db_files_to_remove = \
         build_blast_db_from_fasta_path(db_path)
    created_temp_paths += db_files_to_remove
    
    # split the fasta files and get the list of resulting files        
    tmp_fasta_fps =\
     split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix,output_dir)
    created_temp_paths += tmp_fasta_fps
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)

    # generate the list of commands to be pushed out to nodes    
    commands, job_result_filepaths = get_commands(tmp_fasta_fps,db_path,\
     blastall_fp,blastmat_dir,e_value,word_size,num_hits,output_dir,working_dir,\
     command_prefix=None,command_suffix=None)
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
        write_merge_map_file_blast(job_result_filepaths,output_dir,\
            merge_map_filepath,input_file_basename)
        created_temp_paths.append(merge_map_filepath)
        
        # Create the filepath listing the temporary files to be deleted,
        # but don't write it yet
        deletion_list_filepath = '%s/deletion_list.txt' % working_dir
        created_temp_paths.append(deletion_list_filepath)
        
        if not poll_directly:
            # Generate the command to run the poller, and the list of temp files
            # created by the poller
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
        submit_jobs(cluster_jobs_fp,jobs_fp,job_prefix)
        
    if poll_directly:
        try:
            check_call(poller_command.split())
        except CalledProcessError, e:
            print '**Error occuring when calling the poller directly. '+\
            'Jobs may have been submitted, but are not being polled.'
            print str(e)
            exit(-1)