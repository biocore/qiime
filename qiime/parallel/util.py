#!/usr/bin/env python
#util.py

from __future__ import division
from random import choice
from os import popen, system, getenv, mkdir
from subprocess import Popen, PIPE, STDOUT
from os.path import split
from cogent.parse.fasta import MinimalFastaParser

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"] 
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

RANDOM_JOB_PREFIX_CHARS = "abcdefghigklmnopqrstuvwxyz"
RANDOM_JOB_PREFIX_CHARS += RANDOM_JOB_PREFIX_CHARS.upper()
RANDOM_JOB_PREFIX_CHARS += "0123456790"

def split_fasta(infile, seqs_per_file, outfile_prefix, working_dir=''):
    """ Split infile into files with seqs_per_file sequences in each
    
        infile: list of fasta lines or open file object
        seqs_per_file: the number of sequences to include in each file
        out_fileprefix: string used to create output filepath - output filepaths
         are <out_prefix>.<i>.fasta where i runs from 0 to number of output files
        working_dir: directory to prepend to temp filepaths (defaults to 
         empty string -- files written to cwd)
         
        List of output filepaths is returned.
    
    """
    seq_counter = 0
    out_files = []
    if working_dir and not working_dir.endswith('/'):
        working_dir += '/'
        try:
            mkdir(working_dir)
        except OSError:
            pass
    
    for seq_id,seq in MinimalFastaParser(infile):
        if seq_counter == 0:
            current_out_fp = '%s%s.%d.fasta' \
              % (working_dir,outfile_prefix,len(out_files))
            current_out_file = open(current_out_fp, 'w')
            out_files.append(current_out_fp)
        current_out_file.write('>%s\n%s\n' % (seq_id, seq))
        seq_counter += 1
        
        if seq_counter == seqs_per_file:
            current_out_file.close()
            seq_counter = 0
            
    return out_files

def get_random_job_prefix(fixed_prefix='',max_job_prefix_len=10,\
    leading_trailing_underscores=True):
    """ Return a string to use as job prefix
    
    """

    length = max_job_prefix_len - len(fixed_prefix)
    if leading_trailing_underscores:
        length -= 2 
    
    result = [choice(RANDOM_JOB_PREFIX_CHARS) for i in range(length)]
    if leading_trailing_underscores:
        return fixed_prefix + '_' + ''.join(result) + '_'
    else:
        return fixed_prefix + ''.join(result)

def write_jobs_file(commands,job_prefix=None,jobs_fp=None):
    """ Write commands to jobs_fp and return jobs_fp
    """
    
    if jobs_fp:
        jobs_fp = jobs_fp
    elif job_prefix:
        jobs_fp = job_prefix + 'jobs.txt'
    else:
        jobs_fp = 'jobs.txt'
        
    open(jobs_fp,'w').write('\n'.join(commands))
    
    return jobs_fp
    
def submit_jobs(path_to_cluster_jobs, jobs_fp, job_prefix):
    """ Submit the jobs to the queue using cluster_jobs.py
    """
    cmd = '%s -ms %s %s' % (path_to_cluster_jobs, jobs_fp, job_prefix)
    proc = Popen(cmd,shell=True,universal_newlines=True,\
                 stdout=PIPE,stderr=STDOUT)
    return_value = proc.wait()
    if return_value != 0:
        msg = "\n\n*** Could not start parallel jobs. \n" +\
         "Command run was:\n %s\n" % cmd +\
         "Command returned exit status: %d\n" % return_value +\
         "Stdout/stderr:\n%s\n" % proc.stdout.read()
        raise RuntimeError, msg

def compute_seqs_per_file(input_fasta_fp,num_jobs_to_start):
    """ Compute the number of sequences to include in each split file
    """
    # count the number of sequences in the fasta file
    num_input_seqs = \
     int(popen("egrep -c '^>' %s" % input_fasta_fp).read().strip())
     
    # divide the number of sequences by the number of jobs to start
    result = num_input_seqs/num_jobs_to_start
    
    # if we don't have a perfect split, round up
    if result % 1 != 0:
        result += 1
        
    # return the result as an integer
    return int(result)
    
def build_filepaths_from_filepaths(filepaths,prefix='',directory='',\
    suffix='',replacement=('','')):
    """ Returns a modified list of filepaths
    
        Modifications to each filepath, in order, are:
         1- strip path to convert filepath to filename 
             (e.g., mydir/out.txt => out.txt)
         2- replace all occurrences of replacement[0] with replacement[1]
         3- join directory, prefix, filename, suffix
         
        Order of results corresponds to order of input.
    """

    if directory and not directory.endswith('/'):
        directory += '/'

    results = []
    replace_from = replacement[0]
    replace_to = replacement[1]

    for fp in filepaths:
        file_dir, filename = split(fp)
        filename = filename.replace(replace_from,replace_to)
        result = '%s%s%s%s' % (directory,prefix,filename,suffix)
        results.append(result)
    
    return results


def get_poller_command(python_exe_fp,poller_fp,expected_files_filepath,\
    merge_map_filepath,deletion_list_filepath,seconds_to_sleep,\
    command_prefix='/bin/bash; ',command_suffix='; exit'):
    """Generate command to initiate a poller to monitior/process completed runs
    """
    
    result = '%s %s %s -f %s -m %s -d %s -t %d %s' % \
     (command_prefix,
      python_exe_fp,
      poller_fp,
      expected_files_filepath,
      merge_map_filepath,
      deletion_list_filepath,
      seconds_to_sleep,
      command_suffix)
      
    return result, []
        
def get_rename_command(out_filenames,tmp_output_dir,output_dir):
    """Generate commands to move out_filenames from tmp_output_dir to output_dir
    """
    result = ''
    result_filepaths = []
    for fn in out_filenames:
        tmp_result_filepath = '%s/%s' % (tmp_output_dir,fn)
        result_filepath = '%s/%s' % (output_dir,fn)
        result += \
         '; mv %s %s' % (tmp_result_filepath,result_filepath)
        result_filepaths.append(result_filepath)
    return result, result_filepaths
    
        
def write_filepaths_to_file(job_result_filepaths,expected_files_filepath):
    f = open(expected_files_filepath,'w')
    f.write('\n'.join(job_result_filepaths))
    f.close()  
    
    
def write_merge_map_file_align_seqs(job_result_filepaths,output_dir,\
    merge_map_filepath,input_file_basename):
    
    f = open(merge_map_filepath,'w')
    
    out_filepaths = ['%s/%s_aligned.fasta' % (output_dir,input_file_basename),
                     '%s/%s_failures.fasta' % (output_dir,input_file_basename),
                     '%s/%s_log.txt' % (output_dir,input_file_basename)]
    
    aligned_fps = []
    failures_fps = []
    log_fps = []
    
    for fp in job_result_filepaths:
        if fp.endswith('_aligned.fasta'):
            aligned_fps.append(fp)
        elif fp.endswith('_failures.fasta'):
            failures_fps.append(fp)
        else:
            log_fps.append(fp)
    
    for in_files, out_file in\
     zip([aligned_fps,failures_fps,log_fps],out_filepaths):
        f.write('\t'.join(in_files + [out_file]))
        f.write('\n')
    f.close()
  
def write_merge_map_file_pick_otus(job_result_filepaths,output_dir,\
    merge_map_filepath,input_file_basename):
    
    f = open(merge_map_filepath,'w')
    
    out_filepaths = [\
     '%s/%s_otus.txt' % (output_dir,input_file_basename),
     '%s/%s_otus.log' % (output_dir,input_file_basename)]
    
    otus_fps = []
    log_fps = []
    
    for fp in job_result_filepaths:
        if fp.endswith('_otus.txt'):
            otus_fps.append(fp)
        else:
            log_fps.append(fp)
    
    for in_files, out_file in\
     zip([otus_fps,log_fps],out_filepaths):
        f.write('\t'.join(in_files + [out_file]))
        f.write('\n')
    f.close()
    

def write_merge_map_file_assign_taxonomy(job_result_filepaths,output_dir,\
    merge_map_filepath,input_file_basename):
    
    f = open(merge_map_filepath,'w')
    
    out_filepaths = [\
     '%s/%s_tax_assignments.txt' % (output_dir,input_file_basename),
     '%s/%s_tax_assignments.log' % (output_dir,input_file_basename)]
    
    assignment_fps = []
    log_fps = []
    
    for fp in job_result_filepaths:
        if fp.endswith('_tax_assignments.txt'):
            assignment_fps.append(fp)
        else:
            log_fps.append(fp)
    
    for in_files, out_file in\
     zip([assignment_fps,log_fps],out_filepaths):
        f.write('\t'.join(in_files + [out_file]))
        f.write('\n')
    f.close()
    
def write_merge_map_file_blast(job_result_filepaths,output_dir,\
    merge_map_filepath,input_file_basename):
    
    f = open(merge_map_filepath,'w')
    out_filepath =\
     '%s/%s_blast_out.txt' % (output_dir,input_file_basename)
    f.write('\t'.join(job_result_filepaths + [out_filepath]))
    f.write('\n')
    f.close()
