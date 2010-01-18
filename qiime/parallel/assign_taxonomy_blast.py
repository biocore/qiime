#!/usr/bin/env python
# Author: Greg Caporaso (gregcaporaso@gmail.com)
# parallel_assign_taxonomy_blast.py

""" Description
File created on 13 Sep 2009.

"""
from __future__ import division
from optparse import OptionParser
from os import popen, system, mkdir, makedirs
from os.path import split, splitext
from subprocess import check_call, CalledProcessError
from cogent.app.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_poller_command, get_rename_command, write_filepaths_to_file,\
    write_merge_map_file_assign_taxonomy
from qiime.util import qiime_config

def get_job_commands(python_exe_fp,assign_taxonomy_fp,id_to_taxonomy_fp,\
    e_value,blast_db,\
    blastmat_path,fasta_fps,output_dir,working_dir,\
    command_prefix=None,command_suffix=None):
    """Generate BlastTaxonAssiger classifier commands to be submitted to cluster
    """
    # Create basenames for each of the output files. These will be filled
    # in to create the full list of files created by all of the runs.
    out_filenames = [job_prefix + '.%d_tax_assignments.log', 
                     job_prefix + '.%d_tax_assignments.txt']

    command_prefix = command_prefix or\
     '/bin/bash; cd %s; export BLASTMAT=%s;' \
       % (working_dir,blastmat_path)
    command_suffix = command_suffix or\
     '; exit'
    
    commands = []
    result_filepaths = []
    
    for i,fasta_fp in enumerate(fasta_fps):
        # Each run ends with moving the output file from the tmp dir to
        # the output_dir. Build the command to perform the move here.
        rename_command, current_result_filepaths = get_rename_command(\
         [fn % i for fn in out_filenames],working_dir,output_dir)
        result_filepaths += current_result_filepaths
        
        command = '%s %s %s -o %s -m blast -e %s -b %s -i %s -t %s %s %s' %\
         (command_prefix,\
          python_exe_fp,\
          assign_taxonomy_fp,\
          working_dir,
          e_value,
          blast_db,
          fasta_fp,
          id_to_taxonomy_fp,
          rename_command,
          command_suffix)
          
        commands.append(command)
        
    return commands, result_filepaths

usage_str = """usage: %prog [options] {-i INPUT_FP -o OUTPUT_DIR -t ID_TO_TAXONOMY_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
 BE SURE TO SPECIFY FULL PATHS!
 python Qiime/qiime/parallel/assign_taxonomy_blast.py -O 5 -i /home/caporaso/at_inseqs.fasta -t /home/caporaso/at_id_to_taxonomy.txt -r /home/caporaso/at_refseqs.fasta -o /home/caporaso/out/
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
    
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='full path to store output files '+\
           '[REQUIRED]')

    parser.add_option('-r','--reference_seqs_fp',action='store',\
        help='Ref seqs to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]')
           
    parser.add_option('-b', '--blast_db',
        help='Database to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]')

    parser.add_option('-e', '--e_value', type='float',
        help='Maximum e-value to record an assignment, only used for blast '
        'method [default: %default]',default=0.001)
           
    parser.add_option('-B','--blastmat_dir',action='store',\
           type='string',help='full path to directory containing '+\
           'blastmat file [default: %default]',\
           default=qiime_config['blastmat_dir'])
           
    # Define parallel-script-specific parameters
    parser.add_option('-N','--assign_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/assign_taxonomy.py [default: %default]',\
           default=qiime_config['assign_taxonomy_fp'])
        
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
           
    parser.add_option('-X','--job_prefix',action='store',\
           type='string',help='job prefix '+\
           '[default: BTA_ + 5 random chars]')
           
    parser.add_option('-Y','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default=qiime_config['python_exe_fp'])
        
    parser.add_option('-Z','--seconds_to_sleep',type='int',\
            help='Number of seconds to sleep between checks for run '+\
            ' completion when polling runs [default: %default]',default=60)
                             
    opts,args = parser.parse_args()
    
    required_options = ['input_fasta_fp','id_to_taxonomy_fp','output_dir']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
            
    if not (opts.reference_seqs_fp or opts.blast_db):
        parser.error('Either a blast db (via -b) or a collection of '
                     'reference sequences (via -r) must be passed to '
                     'assign taxonomy using blast.')

    return opts,args
        
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
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
     e_value,blast_db,\
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