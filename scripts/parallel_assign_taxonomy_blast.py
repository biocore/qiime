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
 

from qiime.util import make_option
from os import popen, system, mkdir, makedirs
from os.path import split, splitext, join
from subprocess import check_call, CalledProcessError
from qiime.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_poller_command, get_rename_command, write_filepaths_to_file,\
    write_merge_map_file_assign_taxonomy
from qiime.parallel.assign_taxonomy_blast import get_job_commands
from qiime.util import parse_command_line_parameters, get_options_lookup,\
    get_qiime_scripts_dir, load_qiime_config

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info={}

script_info['brief_description']="""Parallel taxonomy assignment using BLAST"""
script_info['script_description']="""This script performs like the assign_taxonomy.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Example""","""Assign taxonomy to all sequences in the input file (-i) via five (-O) independent jobs using BLAST with the id to taxonomy mapping file (-t) and reference sequence template file (-r), and write the results (-o) to /home/qiime_user/out/. BE SURE TO SPECIFY FULL PATHS!""","""%prog -O 5 -i /home/qiime_user/inseqs.fasta -t /home/qiime_user/at_id_to_taxonomy.txt -r /home/qiime_user/at_refseqs.fasta -o /home/qiime_user/out/"""))

script_info['output_description']="""Mapping of sequence identifiers to taxonomy and quality scores."""

script_info['required_options'] = [\
 make_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp [REQUIRED]'),\
 make_option('-t','--id_to_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'id_to_taxonomy mapping file [REQUIRED]'),\
 make_option('-o','--output_dir',action='store',\
           type='string',help='full path to store output files '+\
           '[REQUIRED]')
]

script_info['optional_options'] = [\
 make_option('-r','--reference_seqs_fp',action='store',\
        help='Ref seqs to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]'),\
 make_option('-b', '--blast_db',
        help='Database to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]'),\
 make_option('-e', '--e_value', type='float',
        help='Maximum e-value to record an assignment, only used for blast '
        'method [default: %default]',default=0.001),\
 make_option('-B','--blastmat_dir',action='store',\
           type='string',help='full path to directory containing '+\
           'blastmat file [default: %default]',\
           default=qiime_config['blastmat_dir']),\
 make_option('-N','--assign_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'scripts/assign_taxonomy.py [default: %default]',\
           default=join(get_qiime_scripts_dir(),'assign_taxonomy.py')),\
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
            
    if not (opts.reference_seqs_fp or opts.blast_db):
        option_parser.error('Either a blast db (via -b) or a collection of '
                     'reference sequences (via -r) must be passed to '
                     'assign taxonomy using blast.')
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    assign_taxonomy_fp = opts.assign_taxonomy_fp
    path_to_cluster_jobs = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    e_value = opts.e_value
    blast_db = opts.blast_db
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    reference_seqs_fp = opts.reference_seqs_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    blastmat_fp = opts.blastmat_dir
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
    # or a random string beginning with BTA (BlastTaxonAssigner)
    job_prefix = opts.job_prefix or get_random_job_prefix('BTA')
    
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
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # Build the blast database from the reference_seqs_fp -- all procs
    # will then access one db rather than create one per proc
    if not blast_db:
        blast_db, db_files_to_remove = \
             build_blast_db_from_fasta_path(reference_seqs_fp)
        created_temp_paths += db_files_to_remove
    else:
        db_files_to_remove = []
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix,output_dir)
    created_temp_paths += tmp_fasta_fps
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)
    
    # generate the list of commands to be pushed out to nodes
    commands, job_result_filepaths  = \
     get_job_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,\
     e_value,blast_db,job_prefix,\
     blastmat_fp,tmp_fasta_fps,output_dir,working_dir)
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
        write_merge_map_file_assign_taxonomy(job_result_filepaths,output_dir,\
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
        submit_jobs(path_to_cluster_jobs,jobs_fp,job_prefix)
    
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