#!/usr/bin/env python
# File created on 14 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.parallel.util import ParallelWrapper


class ParallelMultipleRarefactions(ParallelWrapper):
    _script_name = "single_rarefaction.py"
    _job_prefix = 'RARIF'
    _input_splitter = ParallelWrapper._input_existing_filepaths

    def _identify_files_to_remove(self, job_result_filepaths, params):
        """ The output of the individual jobs are the files we want to keep
        """
        return []

    def _get_job_commands(self,
                          input_fp,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate rarefaction diversity commands to be submitted to cluster
        """
        # Create data for each run (depth, output_fn)
        min_seqs = params['min']
        max_seqs = params['max']
        step = params['step']
        num_reps = params['num_reps']
        run_parameters = []
        for num_seqs in range(min_seqs, max_seqs + 1, step):
            for rep_num in range(num_reps):
                run_parameters.append((
                    num_seqs, 'rarefaction_%d_%d.biom' % (num_seqs, rep_num)))

        commands = []
        result_filepaths = []

        if params['suppress_lineages_included']:
            lineages_included_str = '--suppress_lineages_included'
        else:
            lineages_included_str = ''

        if params['subsample_multinomial']:
            subsample_multinomial_str = '--subsample_multinomial'
        else:
            subsample_multinomial_str = ''

        for depth, output_fn in run_parameters:
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths =\
                self._get_rename_command([output_fn], working_dir, output_dir)
            result_filepaths += current_result_filepaths

            command = '%s %s -i %s -o %s %s %s -d %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 input_fp,
                 working_dir + '/' + output_fn,
                 lineages_included_str,
                 subsample_multinomial_str,
                 depth,
                 rename_command,
                 command_suffix)

            commands.append(command)

        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)
        return commands, result_filepaths
