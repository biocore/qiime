#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import join, split
from qiime.parallel.util import ParallelWrapper


class ParallelAlphaDiversity(ParallelWrapper):
    _script_name = "alpha_diversity.py"
    _job_prefix = 'ALDIV'
    _input_splitter = ParallelWrapper._input_existing_filepaths

    def _identify_files_to_remove(self, job_result_filepaths, params):
        """ The output of the individual jobs are the files we want to keep
        """
        return []

    def _get_job_commands(self,
                          input_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate alpha diversity commands to be submitted to cluster
        """
        commands = []
        result_filepaths = []

        if params['tree_path']:
            tree_str = '-t %s' % params['tree_path']
        else:
            tree_str = ''

        for input_fp in input_fps:
            input_path, input_fn = split(input_fp)
            output_fn = 'alpha_%s' % input_fn
            output_fn = output_fn.replace('.biom', '.txt')
            temp_fp = join(working_dir, output_fn)
            rename_command, current_result_filepaths =\
                self._get_rename_command([output_fn], working_dir, output_dir)
            result_filepaths += current_result_filepaths

            command = '%s %s -i %s -o %s %s -m %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 input_fp,
                 temp_fp,
                 tree_str,
                 params['metrics'],
                 rename_command,
                 command_suffix)

            commands.append(command)

        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)

        return commands, result_filepaths
