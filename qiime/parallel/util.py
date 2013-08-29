#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "Jai Ram Rideout", 
               "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from math import ceil
from os.path import split, splitext, join
from os import makedirs, mkdir
from random import choice
from cogent.parse.fasta import MinimalFastaParser
from qiime.split import split_fasta
from qiime.util import (load_qiime_config,
                        get_qiime_scripts_dir,
                        qiime_system_call,
                        count_seqs)

qiime_config = load_qiime_config()

RANDOM_JOB_PREFIX_CHARS = "abcdefghigklmnopqrstuvwxyz"
RANDOM_JOB_PREFIX_CHARS += RANDOM_JOB_PREFIX_CHARS.upper()
RANDOM_JOB_PREFIX_CHARS += "0123456790"

class ParallelWrapper(object):
    """
    """
    
    def __init__(self,
                 python_exe_fp=qiime_config['python_exe_fp'],
                 cluster_jobs_fp=qiime_config['cluster_jobs_fp'],
                 jobs_to_start=int(qiime_config['jobs_to_start']),
                 poller_fp=join(get_qiime_scripts_dir(),'poller.py'),
                 retain_temp_files=False,
                 suppress_polling=False,
                 seconds_to_sleep=int(qiime_config['seconds_to_sleep'])):
        """  """
        
        self._python_exe_fp = python_exe_fp
        self._cluster_jobs_fp = cluster_jobs_fp
        self._jobs_to_start = jobs_to_start
        self._poller_fp = poller_fp
        self._retain_temp_files = retain_temp_files
        self._suppress_polling = suppress_polling
        self._seconds_to_sleep = seconds_to_sleep

    def _call_initialization(self,
                             input_fp,
                             output_dir,
                             params,
                             job_prefix,
                             poll_directly,
                             suppress_submit_jobs):
        """ Called as the first step in __call__.
        """
        pass

    def _call_cleanup(self,
                      input_fp,
                      output_dir,
                      params,
                      job_prefix,
                      poll_directly,
                      suppress_submit_jobs):
        """ Called as the last step in __call__.
        """
        pass

    def __call__(self,
                 input_fp,
                 output_dir,
                 params,
                 job_prefix=None,
                 poll_directly=False,
                 suppress_submit_jobs=False):
        """ """
        # Generate a list of files and directories that will be cleaned up
        self.files_to_remove = []
        
        # Allow the user to override the default job_prefix (defined by the 
        # base classes)
        if job_prefix is None:
            job_prefix = self._get_random_job_prefix(self._job_prefix)
        
        # A temporary output directory is created in output_dir named
        # job_prefix. Output files are then moved from the temporary
        # directory to the output directory when they are complete,
        # allowing a poller to detect when runs complete by the presence
        # of their output files.
        working_dir = join(output_dir,job_prefix)
        try:
            makedirs(working_dir)
            self.files_to_remove.append(working_dir)
        except OSError:
            # working dir already exists
            pass
        
        # Perform any method-specific setup. This should prevent the need to 
        # overwrite __call__
        self._call_initialization(input_fp,
                                  output_dir,
                                  params,
                                  job_prefix,
                                  poll_directly,
                                  suppress_submit_jobs)
    
        # split the input filepath into directory and filename, base filename and
        # extension for use in naming other files
        try:
            input_dir, input_fn = split(input_fp)
            input_file_basename, input_ext = splitext(input_fn)
        except AttributeError:
            ## THIS IS AWFUL - SOME OF THE SCRIPTS PASS A LIST, SO THE
            ## PREVIOUS BLOCK WON'T WORK... WHAT DO WE WANT TO DO?
            input_dir, input_fn = split(input_fp[0])
            input_file_basename, input_ext = splitext(input_fn)
        
        # Split the input file into the individual job input files. Add the
        # individual job files to the files_to_remove list
        input_fps, remove_input_on_completion = self._input_splitter(
                                         input_fp,
                                         params,
                                         self._jobs_to_start,
                                         job_prefix,
                                         working_dir)
        if remove_input_on_completion:
            self.files_to_remove += input_fps
        
        # Perform any method-specific setup (e.g., formatting a BLAST database)
        self._precommand_initiation(input_fp,output_dir,working_dir,params)
        
        # Generate the list of commands to be pushed out to workers 
        # and the list of output files generated by each job.
        commands, job_result_filepaths = self._get_job_commands(input_fps,
                                                                output_dir,
                                                                params,
                                                                job_prefix,
                                                                working_dir)
        self.files_to_remove += \
         self._identify_files_to_remove(job_result_filepaths,params)

        # Generate the output clean-up files
        merge_map_filepath, deletion_list_filepath, expected_files_filepath =\
         self._initialize_output_cleanup_files(job_result_filepaths,
                                               output_dir,
                                               working_dir,
                                               input_file_basename,
                                               params)

        # Set up poller apparatus if the user does not suppress polling
        if not self._suppress_polling:
            poller_command = self._initiate_polling(job_result_filepaths,
                                                    working_dir,
                                                    poll_directly,
                                                    merge_map_filepath,
                                                    deletion_list_filepath,
                                                    expected_files_filepath)
        
        # If the poller should be run in the same way as the other commands
        # (rather than by the current process), add it to the list of commands
        if not poll_directly:
            commands.append(poller_command)
     
        # Build the filepath for the 'jobs script'. Add that file to the 
        # files_to_remove list.
        jobs_fp = join(working_dir,job_prefix + 'jobs.txt')
        self._write_jobs_file(commands,jobs_fp)
        self.files_to_remove.append(jobs_fp)
    
        # submit the jobs file using cluster_jobs, if not suppressed by the
        # user
        if not suppress_submit_jobs:
            stdout, stderr, return_value = self._submit_jobs(
             jobs_fp=jobs_fp, job_prefix=job_prefix)
        
        # If the poller is going to be run by the current process, 
        # start polling
        if poll_directly:
            # IMPORTANT: the following line MUST use qiime_system_call()
            # instead of subprocess.call, .check_call, or .check_output in case
            # we are invoked in a child process with PIPEs (a deadlock will
            # occur otherwise). This can happen if this code is tested by
            # all_tests.py, for example.
            stdout, stderr, return_value = qiime_system_call(poller_command)
            if return_value != 0:
                print '**Error occuring when calling the poller directly. '+\
                'Jobs may have been submitted, but are not being polled.'
                print stderr
                print poller_command
                exit(-1)
        self.files_to_remove = []

        # Perform any method-specific cleanup. This should prevent the need to 
        # overwrite __call__
        self._call_cleanup(input_fp,
                           output_dir,
                           params,
                           job_prefix,
                           poll_directly,
                           suppress_submit_jobs)

    def _initialize_output_cleanup_files(self,
                                         job_result_filepaths,
                                         output_dir,
                                         working_dir,
                                         input_file_basename,
                                         params):
        # Write the mapping file which described how the output files from
        # each job should be merged into the final output files
        merge_map_filepath = '%s/merge_map.txt' % working_dir
        self._write_merge_map_file(input_file_basename,
                                   job_result_filepaths,
                                   params,
                                   output_dir,
                                   merge_map_filepath)
        self.files_to_remove.append(merge_map_filepath)

        # Create the filepath listing the temporary files to be deleted,
        # but don't write it yet
        deletion_list_filepath = '%s/deletion_list.txt' % working_dir
        self.files_to_remove.append(deletion_list_filepath)
        
        # Write the list of files which must exist for the jobs to be 
        # considered complete
        expected_files_filepath = '%s/expected_out_files.txt' % working_dir
        self._write_filepaths_to_file(job_result_filepaths,
                                      expected_files_filepath)
        self.files_to_remove.append(expected_files_filepath)
        
        return (merge_map_filepath,
                deletion_list_filepath,
                expected_files_filepath)

    def _submit_jobs(self,
                     jobs_fp,
                     job_prefix):
        """ Submit the jobs to the queue using cluster_jobs.py
        """
        cmd = '%s -ms %s %s' % (self._cluster_jobs_fp,
                                jobs_fp, 
                                job_prefix)
        stdout, stderr, return_value = qiime_system_call(cmd)
        if return_value != 0:
            msg = "\n\n*** Could not start parallel jobs. \n" +\
             "Command run was:\n %s\n" % cmd +\
             "Command returned exit status: %d\n" % return_value +\
             "Stdout:\n%s\nStderr\n%s\n" % (stdout,stderr)
            raise RuntimeError, msg
        
        # Leave this comments in as they're useful for debugging.
        # print 'Return value: %d\n' % return_value
        # print 'STDOUT: %s\n' % stdout
        # print 'STDERR: %s\n' % stderr
        
        return stdout, stderr, return_value

    def _identify_files_to_remove(self,job_result_filepaths,params):
        """ Select the files to remove: by default remove all files
        """
        return job_result_filepaths
    
    def _precommand_initiation(self,input_fp,output_dir,working_dir,params):
        pass
    
    def _write_merge_map_file(self,
                              input_file_basename,
                              job_result_filepaths,
                              params,
                              output_dir,
                              merge_map_filepath):
        """ Create an empty file by default. Most subclasses will overwrite this 
        """
        open(merge_map_filepath,'w').close()

    def _get_random_job_prefix(self, 
                               fixed_prefix='',
                               max_job_prefix_len=10,\
                               leading_trailing_underscores=True):
        """ Return a string to use as job prefix """

        length = max_job_prefix_len - len(fixed_prefix)
        if leading_trailing_underscores:
            length -= 2 
    
        result = [choice(RANDOM_JOB_PREFIX_CHARS) for i in range(length)]
        if leading_trailing_underscores:
            return fixed_prefix + '_' + ''.join(result) + '_'
        else:
            return fixed_prefix + ''.join(result)

    def _get_job_commands(self,
                          input_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        raise NotImplementedError, "Subclass must override _get_jobs_commands"
    
    def _initiate_polling(self,
                          job_result_filepaths,
                          working_dir,
                          poll_directly,
                          merge_map_filepath,
                          deletion_list_filepath,
                          expected_files_filepath):
        # Generate the command to run the poller, and the list of temp files
        # created by the poller
        if not poll_directly:
            poller_command, poller_result_filepaths =\
             self._get_poller_command(expected_files_filepath,
                                      merge_map_filepath,
                                      deletion_list_filepath)
        else:
            # this 'else' exists only because of the command_prefix/suffix bullshit,
            # which needs to be refactored (related to trac #109)
            poller_command, poller_result_filepaths =\
             self._get_poller_command(expected_files_filepath,
                                      merge_map_filepath,
                                      deletion_list_filepath,
                                      command_prefix='',
                                      command_suffix='')
        
        self.files_to_remove += poller_result_filepaths
        
        if not self._retain_temp_files:
            # If the user wants temp files deleted, now write the list of 
            # temp files to be deleted
            self._write_filepaths_to_file(self.files_to_remove,
                                    deletion_list_filepath)
        else:
            # Otherwise just write an empty file
            self._write_filepaths_to_file([],
                                    deletion_list_filepath)
        
        return poller_command

    def _get_poller_command(self,
                            expected_files_filepath,
                            merge_map_filepath,
                            deletion_list_filepath,
                            command_prefix='/bin/bash; ',
                            command_suffix='; exit'):
        """Generate command to initiate a poller to monitior/process completed runs
        """
    
        result = '%s poller.py -f %s -m %s -d %s -t %d %s' % \
         (command_prefix,
          expected_files_filepath,
          merge_map_filepath,
          deletion_list_filepath,
          self._seconds_to_sleep,
          command_suffix)
      
        return result, []

    def _write_filepaths_to_file(self,
                                 job_result_filepaths,
                                 expected_files_filepath):
        f = open(expected_files_filepath,'w')
        f.write('\n'.join(job_result_filepaths))
        f.close()  
    
    def _get_rename_command(self,
                            out_filenames,
                            tmp_output_dir,
                            output_dir):
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
    
    def _write_jobs_file(self,
                         commands,
                         jobs_fp):
        """ Write commands to jobs_fp and return jobs_fp """
        open(jobs_fp,'w').write('\n'.join(commands))
        return jobs_fp

    def _merge_to_n_commands(self,
                             commands,
                             n,
                             delimiter=' ; ',
                             command_prefix=None,
                             command_suffix=None):
        """ merge a list of commands into n commands 
            
            This is used by parallel wrappers such as alpha_diversity and 
             beta_diversity which perform an operation on a collection of
             input files (opposed to the scripts that split an input file
             into the user-specified number of jobs).
        
        """
        if n < 1:
            raise ValueError, "number of commands (n) must be an integer >= 1"
        
        commands_to_filter = []
        if command_prefix == None:
            command_prefix = '/bin/bash ;'
            commands_to_filter.append('/bin/bash')
        else:
            commands_to_filter += [c.strip() for c in command_prefix.split(';') if c.strip()]
        
        if command_suffix == None:
            command_suffix = '; exit'
            commands_to_filter.append('exit')
        else:
            commands_to_filter += [c.strip() for c in command_suffix.split(';') if c.strip()]
        
        result = []
        commands_per_merged_command = int(ceil((len(commands)/n)))
        # begin iterating through the commands
        cmd_counter = 0
        current_cmds = []
        for command in commands:
            subcommands = [c.strip() for c in command.split(';')]
            current_cmds.append(delimiter.join([s for s in subcommands if s not in commands_to_filter]))
            cmd_counter += 1
        
            if cmd_counter == commands_per_merged_command:
                result.append(delimiter.join(current_cmds))
                current_cmds = []
                cmd_counter = 0
            
        if current_cmds:
            result[-1] = delimiter.join([result[-1]] + current_cmds)
    
        for i,r in enumerate(result):
            r = '%s %s %s' % (command_prefix, r, command_suffix)
            result[i] = r.strip()
    
        return result
    
    def _compute_seqs_per_file(self,
                               input_fasta_fp,
                               num_jobs_to_start):
        """ Compute the number of sequences to include in each split file
        """
        # count the number of sequences in the fasta file
        num_input_seqs = count_seqs(input_fasta_fp)[0]
     
        # divide the number of sequences by the number of jobs to start
        result = num_input_seqs/num_jobs_to_start
    
        # if we don't have a perfect split, round up
        if result % 1 != 0:
            result += 1
        
        # return the result as an integer
        return int(result)
    
    ####
    # General purpose _input_splitter functions
    ####
    def _split_fasta(self,
                     input_fp,
                     params,
                     jobs_to_start,
                     job_prefix,
                     output_dir):
        # compute the number of sequences that should be included in
        # each file after splitting the input fasta file   
        num_seqs_per_file = self._compute_seqs_per_file(input_fp,jobs_to_start)
     
        # split the fasta files and get the list of resulting files
        tmp_fasta_fps =\
          split_fasta(open(input_fp),num_seqs_per_file,\
          job_prefix,working_dir=output_dir)
        
        return tmp_fasta_fps, True
    
    def _input_existing_filepaths(self,
                                  input_fps,
                                  params,
                                  jobs_to_start,
                                  job_prefix,
                                  output_dir):
        return input_fps, False

class BufferedWriter():
    """A file like object that delays writing to file without keeping an open filehandle

    This class comes useful in scenarios were potentially many open fhs are needed
    (e.g. during splitting of inputs for parallelization). Since
    each OS limits the max number of open fh at any time, we provide a fh like class that
    can be used much like a regular (writable) fh, but without keeping the fh open permanently.
    Using a larger buffer size speeds up the program by using less of the expensive open/close 
    IO operations.    
    """

    def __init__(self, filename, buf_size=100):
        """
        filename: name of file to write to in append mode 
        
        buf_size: buffer size in chunks. Each write operations counts as one chunk.
        """
    
        if(buf_size<1):
            raise ValueError("Invalid buf_size. Must be 1 or larger.")

        self.buffer = []
        self.buf_size = buf_size
        self.filename = filename

        #touch the file
        fh = open(self.filename, "w")
        fh.close()

    def __del__(self):
        self._flush()

    def close(self):
        self._flush()
    
    def write(self, line):
        """write line to BufferedWriter"""

        self.buffer.append(line)
        if (len(self.buffer) > self.buf_size):
            self._flush()

    def _flush(self):
        """Write buffer to file"""

        fh = open(self.filename, "a")
        fh.write("".join(self.buffer))
        fh.close()

        self.buffer = []
