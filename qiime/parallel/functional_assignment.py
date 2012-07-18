#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.parallel.pick_otus import ParallelPickOtus

class ParallelFunctionAssignerUsearch(ParallelPickOtus):
    _script_name = 'functional_assignment.py'
    _job_prefix = 'FUAS'

    def _identify_files_to_remove(self,job_result_filepaths,params):
        """ Select the files to remove: by default remove all files
        """
        # save the .uc files
        result =\
             [fp for fp in job_result_filepaths if not fp.endswith('.uc')]
        return result
    
    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_fmap.txt',
                         job_prefix + '.%s_failures.txt',
                         job_prefix + '.%s.uc',
                         job_prefix + '.%s.bl6']
    
        # Create lists to store the results
        commands = []
        result_filepaths = []
        
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
             '%s %s -i %s -r %s -m usearch -o %s --min_percent_id %s --max_accepts %d --max_rejects %d --queryalnfract %f --targetalnfract %f --evalue %e %s %s' %\
             (command_prefix,
              self._script_name,
              fasta_fp,
              params['refseqs_fp'],
              working_dir,
              params['min_percent_id'],
              params['max_accepts'],
              params['max_rejects'],
              params['queryalnfract'],
              params['targetalnfract'],
              params['evalue'],
              rename_command,
              command_suffix)

            commands.append(command)
        return commands, result_filepaths

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
        blast6_fps = []
    
        out_filepaths = [
         '%s/%s_fmap.txt' % (output_dir,input_file_basename),
         '%s/%s_fmap.log' % (output_dir,input_file_basename),
         '%s/%s_failures.txt' % (output_dir,input_file_basename),
         '%s/%s.bl6' % (output_dir,input_file_basename)]
        in_filepaths = [otus_fps,log_fps,failures_fps,blast6_fps]
    
        for fp in job_result_filepaths:
            if fp.endswith('_fmap.txt'):
                otus_fps.append(fp)
            elif fp.endswith('_fmap.log'):
                log_fps.append(fp)
            elif fp.endswith('_failures.txt'):
                failures_fps.append(fp)
            elif fp.endswith('.bl6'):
                blast6_fps.append(fp)
            else:
                pass
    
        for in_files, out_file in\
         zip(in_filepaths,out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()

