#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"
 

from optparse import make_option
from os import popen, system, mkdir, makedirs
from os.path import split, splitext
from subprocess import check_call, CalledProcessError
from cogent.app.util import get_tmp_filename
from qiime.util import parse_command_line_parameters
from qiime.parallel.assign_taxonomy_rdp import get_commands
from qiime.parallel.util import split_fasta, get_random_job_prefix, write_jobs_file,\
    submit_jobs, compute_seqs_per_file, build_filepaths_from_filepaths,\
    get_poller_command, write_filepaths_to_file,\
    write_merge_map_file_assign_taxonomy
from qiime.util import load_qiime_config

script_description = """ """

script_usage = """Split the input file (-i) into five jobs (-O), start them,
 and write the results (-o) to /home/caporaso/out:

 parallel_assign_taxonomy_rdp.py -O 5 -i 10_seq.fasta -o out"""


qiime_config = load_qiime_config()

required_options = [\
 make_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp [REQUIRED]'),\
 make_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files '+\
           '[REQUIRED]'),\
]

optional_options = [\
 make_option('-r','--rdp_classifier_fp',action='store',\
           type='string',help='full path to rdp classifier jar file '+\
           '[default: %default]',\
           default=qiime_config['rdp_classifier_fp']),\
 make_option('-c','--confidence',action='store',\
          type='float',help='Minimum confidence to'+\
          ' record an assignment [default: %default]',default=0.80),\
 make_option('-N','--assign_taxonomy_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/assign_taxonomy.py [default: %default]',\
           default=qiime_config['assign_taxonomy_fp']),\
 make_option('-O','--jobs_to_start',action='store',type='int',\
            help='Number of jobs to start [default: %default]',default=24),\
 make_option('-P','--poller_fp',action='store',\
           type='string',help='full path to '+\
           'qiime/parallel/poller.py [default: %default]',\
           default=qiime_config['poller_fp']),\
 make_option('-R','--retain_temp_files',action='store_true',\
           help='retain temporary files after runs complete '+\
           '(useful for debugging) [default: %default]',\
           default=False),\
 make_option('-S','--suppress_submit_jobs',action='store_true',\
            help='Only split input and write commands file - don\'t submit '+\
            'jobs [default: %default]',default=False),\
 make_option('-T','--poll_directly',action='store_true',\
            help='Poll directly for job completion rather than running '+\
            'poller as a separate job. If -T is specified this script will '+\
            'not return until all jobs have completed. [default: %default]',\
            default=False),\
 make_option('-U','--cluster_jobs_fp',action='store',\
            type='string',help='path to cluster_jobs.py script ' +\
            ' [default: %default]',\
            default=qiime_config['cluster_jobs_fp']),\
 make_option('-W','--suppress_polling',action='store_true',
           help='suppress polling of jobs and merging of results '+\
           'upon completion [default: %default]',\
           default=False),\
 make_option('-X','--job_prefix',action='store',\
           type='string',help='job prefix '+\
           '[default: RDP_ + 5 random chars]'),\
 make_option('-Y','--python_exe_fp',action='store',\
           type='string',help='full path to python '+\
           'executable [default: %default]',\
           default=qiime_config['python_exe_fp']),\
 make_option('-Z','--seconds_to_sleep',type='int',\
            help='Number of seconds to sleep between checks for run '+\
            ' completion when polling runs [default: %default]',default=60)
]



def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    
    # create local copies of command-line options
    python_exe_fp = opts.python_exe_fp
    assign_taxonomy_fp = opts.assign_taxonomy_fp
    confidence = opts.confidence
    rdp_classifier_fp = opts.rdp_classifier_fp
    cluster_jobs_fp = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    poll_directly = opts.poll_directly

    try:
        makedirs(output_dir)
    except OSError:
        # output dir already exists
        pass

    created_temp_paths = []
    
    # split the input filepath into directory and filename, base filename and
    # extension
    input_dir, input_fasta_fn = split(input_fasta_fp)
    input_file_basename, input_fasta_ext = splitext(input_fasta_fn)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with RDP
    job_prefix = opts.job_prefix or get_random_job_prefix('RDP')

    # A temporary output directory is created in output_dir named
    # job_prefix. Output files are then moved from the temporary 
    # directory to the output directory when they are complete, allowing
    # a poller to detect when runs complete by the presence of their
    # output files.
    working_dir = '%s/%s' % (output_dir,job_prefix)
    try:
        mkdir(working_dir)
        created_temp_paths.append(working_dir)
    except OSError:
        # working dir already exists
        pass
    
    # compute the number of sequences that should be included in
    # each file after splitting the input fasta file   
    num_seqs_per_file = compute_seqs_per_file(input_fasta_fp,jobs_to_start)
     
    # split the fasta files and get the list of resulting files
    tmp_fasta_fps =\
      split_fasta(open(input_fasta_fp),num_seqs_per_file,job_prefix,output_dir)
    created_temp_paths += tmp_fasta_fps
    
    # build the filepath for the 'jobs script'
    jobs_fp = '%s/%sjobs.txt' % (output_dir, job_prefix)
    created_temp_paths.append(jobs_fp)

    # generate the list of commands to be pushed out to nodes
    commands, job_result_filepaths = \
     get_commands(python_exe_fp,assign_taxonomy_fp,confidence,job_prefix,\
     tmp_fasta_fps,rdp_classifier_fp,output_dir,working_dir)
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