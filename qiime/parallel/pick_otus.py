#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jens Reeder"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from cogent.app.formatdb import build_blast_db_from_fasta_path
from qiime.parallel.util import ParallelWrapper
from qiime.parallel.poller import basic_process_run_results_f
from cogent.parse.fasta import MinimalFastaParser
from os.path import basename
from datetime import datetime

class ParallelPickOtus(ParallelWrapper):
    _script_name = "pick_otus.py"
    _job_prefix = 'POTU'
    _input_splitter = ParallelWrapper._split_fasta
    _process_run_results_f =\
         'qiime.parallel.pick_otus.parallel_pick_otus_process_run_results_f'
    
    def _get_poller_command(self,
                            expected_files_filepath,
                            merge_map_filepath,
                            deletion_list_filepath,
                            command_prefix='/bin/bash; ',
                            command_suffix='; exit'):
        """Generate command to initiate a poller to monitior/process completed runs
        """

        result = '%s poller.py -f %s -p %s -m %s -d %s -t %d %s' % \
         (command_prefix,
          expected_files_filepath,
          self._process_run_results_f,
          merge_map_filepath,
          deletion_list_filepath,
          self._seconds_to_sleep,
          command_suffix)
  
        return result, []

    def _write_merge_map_file(self,
                              input_file_basename,
                              job_result_filepaths,
                              params,
                              output_dir,
                              merge_map_filepath):
        """ 
        """
        f = open(merge_map_filepath,'w')
    
        otus_fps = []
        log_fps = []
        failures_fps = []
        
        # determine if the otu picking method generates
        # failures files
        failures = False
        for fp in job_result_filepaths:
            if fp.endswith('_failures.txt'):
                failures = True
                break
        
        if not failures:
            out_filepaths = [
             '%s/%s_otus.txt' % (output_dir,input_file_basename),
             '%s/%s_otus.log' % (output_dir,input_file_basename)]
            in_filepaths = [otus_fps,log_fps]
        else:
            out_filepaths = [
             '%s/%s_otus.txt' % (output_dir,input_file_basename),
             '%s/%s_otus.log' % (output_dir,input_file_basename),
             '%s/%s_failures.txt' % (output_dir,input_file_basename)]
            in_filepaths = [otus_fps,log_fps,failures_fps]
    
        for fp in job_result_filepaths:
            if fp.endswith('_otus.txt'):
                otus_fps.append(fp)
            elif fp.endswith('_otus.log'):
                log_fps.append(fp)
            elif fp.endswith('_failures.txt'):
                failures_fps.append(fp)
            else:
                pass
    
        for in_files, out_file in\
         zip(in_filepaths,out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()

class ParallelPickOtusUclustRef(ParallelPickOtus):
    
    def _identify_files_to_remove(self,job_result_filepaths,params):
        """ Select the files to remove: by default remove all files
        """
        if params['save_uc_files']:
            # keep any .uc files that get created
            result =\
             [fp for fp in job_result_filepaths if not fp.endswith('.uc')]
        else:
            result = [job_result_filepaths]
        
        return result
    
    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be run
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_otus.log', 
                         job_prefix + '.%d_otus.txt',
                         job_prefix + '.%s_failures.txt']
    
        # Create lists to store the results
        commands = []
        result_filepaths = []
    
        if params['enable_rev_strand_match']:
            enable_rev_strand_match_str = '-z'
        else:
            enable_rev_strand_match_str = ''
            
        if params['optimal_uclust']:
            optimal_uclust_str = '-A'
        else:
            optimal_uclust_str = ''
            
        if params['exact_uclust']:
            exact_uclust_str = '-E'
        else:
            exact_uclust_str = ''
            
        if params['stable_sort']:
            stable_sort_str = ''
        else:
            stable_sort_str = '--suppress_uclust_stable_sort'
            
        if params['save_uc_files']:
            save_uc_files_str = ''
            out_filenames += [job_prefix + '%d_clusters.uc']
        else:
            save_uc_files_str = '-d'
        
        # Iterate over the input files
        for i,fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % i for fn in out_filenames],
                working_dir,
                output_dir)
            result_filepaths += current_result_filepaths
            
            command = \
             '%s %s -i %s -r %s -m uclust_ref --suppress_new_clusters -o %s -s %s %s %s %s --max_accepts %s --max_rejects %s --stepwords %d --w %d %s %s %s %s' %\
             (command_prefix,
              self._script_name,\
              fasta_fp,\
              params['refseqs_fp'],\
              working_dir,\
              params['similarity'],\
              enable_rev_strand_match_str,
              optimal_uclust_str,
              exact_uclust_str,
              params['max_accepts'],
              params['max_rejects'],
              params['stepwords'],
              params['word_length'],
              stable_sort_str,
              save_uc_files_str,
              rename_command,
              command_suffix)

            commands.append(command)

        return commands, result_filepaths

class ParallelPickOtusBlast(ParallelPickOtus):

    def _precommand_initiation(self,input_fp,output_dir,working_dir,params):
        if not params['blast_db']:        
            # Build the blast database from the reference_seqs_fp -- all procs
            # will then access one db rather than create one per proc
            blast_db, db_files_to_remove = \
                 build_blast_db_from_fasta_path(params['refseqs_fp'])
            self.files_to_remove += db_files_to_remove
            params['blast_db'] = blast_db

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be submitted to cluster
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_otus.log', 
                         job_prefix + '.%d_otus.txt']
    
        # Create lists to store the results
        commands = []
        result_filepaths = []
    
        # Iterate over the input files
        for i,fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = self._get_rename_command(\
             [fn % i for fn in out_filenames],working_dir,output_dir)
            result_filepaths += current_result_filepaths
            
            command = \
             '%s %s -i %s -b %s -m blast -o %s -e %s -s %s --min_aligned_percent %s %s %s' %\
             (command_prefix,
              self._script_name,
              fasta_fp,
              params['blast_db'],
              working_dir,
              params['max_e_value'],
              params['similarity'],
              params['min_aligned_percent'],
              rename_command,
              command_suffix)
          
            commands.append(command)

        return commands, result_filepaths


def parallel_pick_otus_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles
    
        f: file containing one set of mapping instructions per line
        
        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log
         f1_failures.txt f2_failures.txt f3_failures.txt f_failires.txt
         
        If f contained the two lines above, this function would 
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    lines = list(f)
    # handle catting of log files and failure files
    basic_process_run_results_f([lines[1]])
    try:
        basic_process_run_results_f([lines[2]])
        basic_process_run_results_f([lines[3]])
    except IndexError:
        # no failures files or blast6 were generated (BLAST
        # doesn't create these)
        pass
    # handle merging of otu maps
    fields = lines[0].strip().split()
    infiles_list = fields[:-1]
    out_filepath = fields[-1] 
    try:
        of = open(out_filepath,'w')
    except IOError:
        raise IOError,\
         "Poller can't open final output file: %s" % out_filepath  +\
         "\nLeaving individual jobs output.\n Do you have write access?"

    unique_otu_map = {}
    for fp in infiles_list:
        for line in open(fp):
            fields = line.strip().split()
            try:
                # current otu_id already exists, so append this
                # set of seq_ids
                unique_otu_map[fields[0]] += fields[1:]
            except KeyError:
                # current otu_id has not been seen yet, so 
                # create it with the current set of otus
                unique_otu_map[fields[0]] = fields[1:]
    
    for otu_id, seq_ids in unique_otu_map.items():
        of.write('\t'.join([otu_id] + seq_ids))
        of.write('\n')            
    of.close()
    
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback, 
    # you'll get an error right away rather than going into an infinite loop
    return True


#TODO: 
# - distribut jobs evely according to load in multicore env,
#   if we have a queuing system we could submit all jobs individually
# - ask Greg if there is a more direct way to make the prefix_len available
#   to the input_splitter, as it does not get the params dict. Right now, I 
#   have to overwrite __init__ just to pass this param in.

class ParallelPickOtusTrie(ParallelPickOtus):

    _process_run_results_f =\
        'qiime.parallel.pick_otus.parallel_pick_otus_trie_process_run_results_f'

    def __init__(self,
                 prefix_length=1,
                 *args,
                 **kwargs):
        super(ParallelPickOtusTrie, self).__init__(*args, **kwargs)
        self.prefix_length = prefix_length

        # Maybe use this to distribute across fixed number of jobs with relatively even load 
        self.prefix_counts = {}

    def _split_along_prefix(self, input_fp, jobs_to_start,
                            job_prefix, output_dir):
        """ Split input sequences into sets with identical prefix
        """
        out_files = []
        buffers = {}
        for seq_id, seq in MinimalFastaParser(open(input_fp)):

            if(len(seq) < self.prefix_length):
                raise ValueError("Prefix length must be equal or longer than sequence.\n"
                                 +" Found seq %s with length %d" % (seq_id, len(seq)))
            prefix = seq[:self.prefix_length]

            if (prefix not in buffers):
                # never seen this prefix before
                out_fp = "%s/%s%s" % (output_dir, job_prefix, prefix)
                buffers[prefix] = BufferedWriter(out_fp)
                out_files.append(out_fp)
                self.prefix_counts[prefix] = 0

            self.prefix_counts[prefix] += 1                                
            buffers[prefix].write('>%s\n%s\n' % (seq_id, seq))

        remove_files=True 
        return out_files, remove_files

    _input_splitter = _split_along_prefix

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be run
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = ['%s_otus.log', 
                         '%s_otus.txt']
    
        # Create lists to store the results
        commands = []
        result_filepaths = []
            
        # Iterate over the input files
        for i,fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % basename(fasta_fp) for fn in out_filenames],
                working_dir,
                output_dir)
            result_filepaths += current_result_filepaths
            
            command = \
             '%s %s -i %s -m trie -o %s %s %s' %\
             (command_prefix,
              self._script_name,\
              fasta_fp,\
              working_dir,\
              rename_command,
              command_suffix)

            commands.append(command)

        #TODO: on a queue based cluster system, we should let 
        #      the queue do this and submit all jobs as they are
        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)
        return commands, result_filepaths


#TODO: This class needs tests and maybe should be moved to some of the qiime util modules if generally useful
class BufferedWriter():
    """A file like thing that delays writing to file without keeping an open filehandle

    This class comes useful in scenarios were potentially many open fhs are needed. Since
    each OS limits the max number of open fh at any time, we provide a fh like class that
    can be used much like a regular (writable) fh, but without keeping the fh open permanently.
    """

    def __init__(self, filename, buf_size=100):
        """
        filename: name of file to write to in append mode 
        
        buf_size: buffer size in lines
        """
        self.buffer = []
        self.buf_size = buf_size
        self.filename = filename

        #test open the file
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

def parallel_pick_otus_trie_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles
    
        f: file containing one set of mapping instructions per line
        
        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log
         f1_failures.txt f2_failures.txt f3_failures.txt f_failires.txt
         
        If f contained the two lines above, this function would 
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
         
         Note: this is mostly copied from 
               parallel_pick_otus_process_run_results_f
               The only difference is the way how otu_maps are combined.
               Since each otu_map for the parallel trie otu pickers yields disjoint
               sets, we have to make the otu_ids disjoint as well.
    """
    lines = list(f)
    # handle catting of log files and failure files
    basic_process_run_results_f([lines[1]])
    try:
        basic_process_run_results_f([lines[2]])
        basic_process_run_results_f([lines[3]])
    except IndexError:
        # no failures files or blast6 were generated (BLAST
        # doesn't create these)
        pass
    # handle merging of otu maps
    fields = lines[0].strip().split()
    infiles_list = fields[:-1]
    out_filepath = fields[-1] 
    try:
        of = open(out_filepath,'w')
    except IOError:
        raise IOError,\
         "Poller can't open final output file: %s" % out_filepath  +\
         "\nLeaving individual jobs output.\n Do you have write access?"

    unique_otu_map = {}
    otu_id = 0
    for fp in infiles_list:
        for line in open(fp):
            fields = line.strip().split()
            of.write('\t'.join(["%d" % otu_id] + fields[1:]))
            of.write('\n')            
            otu_id += 1
    of.close()
    
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback, 
    # you'll get an error right away rather than going into an infinite loop
    return True
