#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso","Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from qiime.util import make_option
from qiime.util import parse_command_line_parameters,\
    load_qiime_config, get_options_lookup, get_qiime_scripts_dir
from qiime.parallel.pick_otus_uclust_ref import get_job_commands,\
    get_poller_command
from qiime.parallel.util import split_fasta, get_random_job_prefix,\
    write_jobs_file,submit_jobs, compute_seqs_per_file,\
    build_filepaths_from_filepaths, write_filepaths_to_file,\
    write_merge_map_file_pick_otus
from os import popen, system, makedirs, mkdir
from os.path import split, splitext, join
from subprocess import check_call, CalledProcessError
from qiime.util import get_tmp_filename
from cogent.app.formatdb import build_blast_db_from_fasta_path

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Parallel pick otus using uclust_ref"""
script_info['script_description']="""This script works like the pick_otus.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example""","""Pick OTUs with uclust_ref by searching /home/qiime/inseqs.fasta against /home/qiime/refseqs.fasta and write the output to the /home/qiime/out/ directory.""","""%prog -i /home/qiime/inseqs.fasta -r /home/qiime/refseqs.fasta -o /home/qiime/out/"""))
script_info['output_description']="""The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an OTU identifier which is the reference sequence identifier, and the remaining fields correspond to sequence identifiers assigned to that OTU. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file."""

script_info['required_options'] = [\
    make_option('-i','--input_fasta_fp',action='store',\
           type='string',help='full path to '+\
           'input_fasta_fp'),
           
    make_option('-o','--output_dir',action='store',\
           type='string',help='path to store output files'),
          
    make_option('-r','--refseqs_fp',action='store',\
           type='string',help='full path to '+\
           'reference collection')
]

script_info['optional_options'] = [\
         
    make_option('-s','--similarity',action='store',\
          type='float',help='Sequence similarity '+\
          'threshold [default: %default]',default=0.97),\
    make_option('-z', '--enable_rev_strand_match', action='store_true',
        default=False,
        help=('Enable reverse strand matching for uclust otu picking, '
              'will double the amount of memory used. [default: %default]')),
    make_option('-A','--optimal_uclust', action='store_true', 
              default=False,
              help=('Pass the --optimal flag to uclust for uclust otu'
              ' picking. [default: %default]')),
    make_option('-E','--exact_uclust', action='store_true', 
              default=False,
              help=('Pass the --exact flag to uclust for uclust otu'
              ' picking. [default: %default]')),
    make_option('--max_accepts',type='int',default=20,
              help="max_accepts value to uclust and uclust_ref [default: %default]"),
    make_option('--max_rejects',type='int',default=500,
              help="max_rejects value to uclust and uclust_ref [default: %default]"),
    make_option('--stepwords',type='int',default=20,
              help="stepwords value to uclust and uclust_ref [default: %default]"),
    make_option('--word_length',type='int',default=12,
              help="w value to uclust and uclust_ref [default: %default]"),
    make_option('--uclust_stable_sort',default=True,action='store_true',
              help=("Deprecated: stable sort enabled by default, pass "
                  "--uclust_suppress_stable_sort to disable [default: %default]")),         
    make_option('--suppress_uclust_stable_sort',default=False,action='store_true',
                  help=("Don't pass --stable-sort to uclust [default: %default]")),
    make_option('-d', '--save_uc_files', default=True, action='store_false',
              help=("Enable preservation of intermediate uclust (.uc) files "
              "that are used to generate clusters via uclust. "
              "[default: %default]")),
    #Define parallel-script-specific parameters
    make_option('-N','--pick_otus_fp',action='store',\
           type='string',help='full path to '+\
           'scripts/pick_otus.py [default: %default]',\
           default=join(get_qiime_scripts_dir(),'pick_otus.py')),\
        
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
    python_exe_fp = opts.python_exe_fp
    pick_otus_fp = opts.pick_otus_fp
    refseqs_fp = opts.refseqs_fp
    cluster_jobs_fp = opts.cluster_jobs_fp
    input_fasta_fp = opts.input_fasta_fp 
    jobs_to_start = opts.jobs_to_start
    output_dir = opts.output_dir
    poller_fp = opts.poller_fp
    retain_temp_files = opts.retain_temp_files
    suppress_polling = opts.suppress_polling
    seconds_to_sleep = opts.seconds_to_sleep
    similarity = opts.similarity
    poll_directly = opts.poll_directly
    uclust_stable_sort = not opts.suppress_uclust_stable_sort
    save_uc_files = opts.save_uc_files
    
    enable_rev_strand_match = opts.enable_rev_strand_match
    optimal_uclust = opts.optimal_uclust
    exact_uclust = opts.exact_uclust
    max_accepts = opts.max_accepts
    max_rejects = opts.max_rejects
    stepwords = opts.stepwords
    word_length = opts.word_length

    created_temp_paths = []
    
    # split the input filepath into directory and filename, base filename and
    # extension
    input_dir, input_fasta_fn = split(input_fasta_fp)
    input_file_basename, input_fasta_ext = splitext(input_fasta_fn)
    
    # set the job_prefix either based on what the user passed in,
    # or a random string beginning with RDP
    job_prefix = opts.job_prefix or get_random_job_prefix('POTU')
    
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
     get_job_commands(python_exe_fp,pick_otus_fp,tmp_fasta_fps,
     output_dir,refseqs_fp,job_prefix,working_dir,similarity,
     enable_rev_strand_match,optimal_uclust,exact_uclust,max_accepts,max_rejects,
     stepwords, word_length, uclust_stable_sort, save_uc_files)
    if save_uc_files:
        # keep any .uc files that get created
        created_temp_paths +=\
         [fp for fp in job_result_filepaths if not fp.endswith('.uc')]
    else:
        created_temp_paths += [job_result_filepaths]

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
        process_run_results_f =\
         'qiime.parallel.pick_otus_uclust_ref.parallel_uclust_ref_process_run_results_f'
        write_merge_map_file_pick_otus(job_result_filepaths,output_dir,\
            merge_map_filepath,input_file_basename,failures=True)
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
             merge_map_filepath,deletion_list_filepath,process_run_results_f,\
             seconds_to_sleep=seconds_to_sleep)
            created_temp_paths += poller_result_filepaths
            # append the poller command to the list of job commands
            commands.append(poller_command)
        else:
            poller_command, poller_result_filepaths =\
             get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
             merge_map_filepath,deletion_list_filepath,process_run_results_f,\
             seconds_to_sleep=seconds_to_sleep,command_prefix='',command_suffix='')
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
