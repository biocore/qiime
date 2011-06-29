#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from os import makedirs, system
from os.path import exists, split, splitext
from subprocess import check_call, CalledProcessError
from optparse import OptionParser
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_poller_command, write_filepaths_to_file,\
    write_merge_map_file_blast
from qiime.parallel.blast import get_commands
from qiime.util import load_qiime_config, get_options_lookup

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Parallel BLAST"""
script_info['script_description']="""This script for performing blast while making use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example""","""BLAST /home/qiime_user/10_seq.fasta (-i) via three (-O) independent jobs against a blast database created from /home/qiime_user/1000_seq.fasta (-r). Store the results in /home/qiime_user/bla_out/ (-o).""","""%prog -i /home/qiime_user/10_seq.fasta -r /home/qiime_user/1000_seq.fasta -O 3 -o /home/qiime_user/bla_out/"""))
script_info['output_description']=""" """
script_info['required_options'] = [\
 make_option('-i','--infile_path',action='store',\
          type='string',dest='infile_path',
          help='Path of sequences to use as queries [REQUIRED]'),\
 make_option('-r','--refseqs_path',action='store',\
          type='string',
            help='Path to fasta sequences to search against' +\
            ' [REQUIRED]'),\
 make_option('-o', '--output_dir', \
        help='name of output directory for blast jobs [REQUIRED]')
]
script_info['optional_options'] = [\
 make_option('-c','--disable_low_complexity_filter',
        default=False,action='store_true',
        help='disable filtering of low-complexity sequences '
             '(i.e., -F F is passed to blast) [default: %default]'),\
 make_option('-e','--e_value',action='store',\
        type='float', default=1e-30, dest='e_value',
        help='E-value threshold for blasts [default: %default]'),\
 make_option('-n','--num_hits',action='store',\
        type='int', default=1, dest='num_hits',
        help='number of hits per query for blast results [default: %default]'),\
 make_option('-w','--word_size',action='store',\
        type='int', default=30, dest='word_size',
        help='word size for blast searches [default: %default]'),\
 make_option('-D', '--suppress_format_blastdb', action='store_true',\
        default=False,help='supress format of blastdb [default: %default]'),\
 make_option('-a','--blastmat_dir',action='store',\
           type='string',help='full path to directory containing '+\
           'blastmat file [default: %default]',\
           default=qiime_config['blastmat_dir']),\
 make_option('-b','--blastall_fp',
        default=qiime_config['blastall_fp'],
        help='Path to blastall [default: %default]'),\
 options_lookup['jobs_to_start'],\
 options_lookup['poller_fp'],\
 options_lookup['retain_temp_files'],\
 options_lookup['suppress_submit_jobs'],\
 options_lookup['poll_directly'],\
 options_lookup['cluster_jobs_fp'],\
 options_lookup['suppress_polling'],\
 options_lookup['job_prefix'],\
 options_lookup['python_exe_fp'],\
 options_lookup['seconds_to_sleep']\
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)    
    # create local copies of command-line options
    input_fasta_fp = opts.infile_path 
    refseqs_path = opts.refseqs_path
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
    disable_low_complexity_filter = opts.disable_low_complexity_filter

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
         build_blast_db_from_fasta_path(refseqs_path)
    created_temp_paths += db_files_to_remove
    
    # split the fasta files and get the list of resulting files        
    tmp_fasta_fps =\
     split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix,output_dir)
    created_temp_paths += tmp_fasta_fps
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)

    # generate the list of commands to be pushed out to nodes    
    commands, job_result_filepaths = get_commands(tmp_fasta_fps,refseqs_path,
     blastall_fp,blastmat_dir,e_value,word_size,num_hits,output_dir,working_dir,
     disable_low_complexity_filter=disable_low_complexity_filter,
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


if __name__ == "__main__":
    main()