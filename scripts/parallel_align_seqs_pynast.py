#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
from qiime.util import make_option
from os.path import split, splitext, join
from os import popen, system, makedirs, mkdir
from subprocess import check_call, CalledProcessError
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.util import load_qiime_config, parse_command_line_parameters,\
    get_options_lookup, get_qiime_scripts_dir
from qiime.align_seqs import compute_min_alignment_length
from pynast.util import pairwise_alignment_methods
from qiime.parallel.util import split_fasta, get_random_job_prefix,\
    write_jobs_file, submit_jobs, compute_seqs_per_file,\
    build_filepaths_from_filepaths, get_poller_command,\
    write_filepaths_to_file, write_merge_map_file_align_seqs
from qiime.parallel.align_seqs_pynast import get_job_commands

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Parallel sequence alignment using PyNAST"""
script_info['script_description']="""A wrapper for the align_seqs.py PyNAST option, intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example""","""Align the input file (-i) against using PyNAST and write the output (-o) to $PWD/pynast_aligned_seqs/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""","""%prog -i $PWD/inseqs.fasta -o $PWD/pynast_aligned_seqs/"""))
script_info['output_description']="""This results in a multiple sequence alignment (FASTA-formatted)."""

script_info['required_options'] = [\
 options_lookup['fasta_as_primary_input'],\
 options_lookup['output_dir']
]

pairwise_alignment_method_choices = pairwise_alignment_methods.keys()
blast_db_default_help =\
 qiime_config['pynast_template_alignment_blastdb'] or \
 'created on-the-fly from template_alignment'

script_info['optional_options'] = [\
 make_option('-a','--pairwise_alignment_method',\
          type='choice',help='Method to use for pairwise alignments'+\
          ' [default: %default]',\
          default='uclust',choices=pairwise_alignment_method_choices),\
 make_option('-d','--blast_db',\
          dest='blast_db',help='Database to blast against'+\
          ' [default: %s]' % blast_db_default_help,
          default=qiime_config['pynast_template_alignment_blastdb']),\
    make_option('-e','--min_length',\
          type='int',help='Minimum sequence '+\
          'length to include in alignment [default: 75% of the'+\
          ' median input sequence length]',\
           default=-1),
 make_option('-p','--min_percent_id',action='store',\
          type='float',help='Minimum percent '+\
          'sequence identity to closest blast hit to include sequence in'+\
          ' alignment [default: %default]',default=75.0),\
 make_option('-N','--align_seqs_fp',\
           type='existing_filepath',help='full path to '+\
           'Qiime/scripts/align_seqs.py [default: %default]',\
           default=join(get_qiime_scripts_dir(),'align_seqs.py')),\
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

# pynast_template_alignment_fp is required only if it is not 
# provided in qiime_config
if qiime_config['pynast_template_alignment_fp']:
    script_info['optional_options'].append(make_option('-t','--template_fp',\
      type='string',dest='template_fp',help='Filepath for '+\
      'template against [default: %default]',
      default=qiime_config['pynast_template_alignment_fp']))
else:
    script_info['required_options'].append(make_option('-t','--template_fp',\
      type='string',dest='template_fp',\
      help='Filepath for template against',
      default=qiime_config['pynast_template_alignment_fp']))

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    align_seqs_fp = opts.align_seqs_fp
    cluster_jobs_fp = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    template_aln_fp = opts.template_fp
    pairwise_alignment_method = opts.pairwise_alignment_method
    min_length = opts.min_length
    min_percent_id = opts.min_percent_id
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    poll_directly = opts.poll_directly

    created_temp_paths = []

    # compute the minimum alignment length if a negative value was
    # provided (the default)
    min_length = opts.min_length
    if min_length < 0:
        min_length = compute_min_alignment_length(open(input_fasta_fp,'U'))

    # split the input filepath into directory and filename, base filename and
    # extension
    input_dir, input_fasta_fn = split(input_fasta_fp)
    input_file_basename, input_fasta_ext = splitext(input_fasta_fn)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with RDP
    job_prefix = opts.job_prefix or get_random_job_prefix('ALIGN')
    
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

    if not opts.blast_db:
        # Build the blast database from the reference_seqs_fp -- all procs
        # will then access one db rather than create one per proc
        blast_db, db_files_to_remove = \
             build_blast_db_from_fasta_path(template_aln_fp,output_dir=working_dir)
        created_temp_paths += db_files_to_remove
    else:
        blast_db = opts.blast_db
    
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
     get_job_commands(python_exe_fp,align_seqs_fp,tmp_fasta_fps,template_aln_fp,\
     pairwise_alignment_method,output_dir,blast_db,\
     min_length,min_percent_id,job_prefix,working_dir)
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
        write_merge_map_file_align_seqs(job_result_filepaths,output_dir,\
            merge_map_filepath,input_file_basename)
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
            created_temp_paths += poller_result_filepaths
            # append the poller command to the list of job commands
            commands.append(poller_command)
        else:
            poller_command, poller_result_filepaths =\
             get_poller_command(python_exe_fp,poller_fp,\
              expected_files_filepath,merge_map_filepath,\
              deletion_list_filepath,seconds_to_sleep=seconds_to_sleep,\
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