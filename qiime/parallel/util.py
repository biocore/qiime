#!/usr/bin/env python
#util.py

from __future__ import division
from random import choice
from os import popen, system
from os.path import split
from cogent.parse.fasta import MinimalFastaParser

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

def write_jobs_file(commands,job_prefix=None,job_fp=None):
    """ Write commands to job_fp and return job_fp
    """
    
    if job_fp:
        job_fp = job_fp
    elif job_prefix:
        job_fp = job_prefix + 'jobs.txt'
    else:
        job_fp = 'jobs.txt'
        
    open(job_fp,'w').write('\n'.join(commands))
    
    return job_fp
    
def submit_jobs(path_to_cluster_jobs, jobs_fp, job_prefix):
    """ Submit the jobs to the queue using cluster_jobs.py
    """
    system('%s -ms %s %s' %\
     (path_to_cluster_jobs, jobs_fp, job_prefix))

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

# if __name__ == '__main__':
#     from sys import argv
#     infilename = argv[1]
#     infile = open(infilename, 'U')
#     seqs_per_file = int(argv[2])
#     split_fasta(infile, seqs_per_file, infilename)
