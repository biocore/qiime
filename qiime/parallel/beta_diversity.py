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

from os.path import join, split, splitext
from skbio.util import create_dir
from biom import load_table

from qiime.parallel.util import ParallelWrapper
from qiime.format import format_distance_matrix


class ParallelBetaDiversity(ParallelWrapper):
    _script_name = "beta_diversity.py"
    _input_splitter = ParallelWrapper._input_existing_filepaths
    _job_prefix = 'BDIV'

    def _identify_files_to_remove(self, job_result_filepaths, params):
        """ The output of the individual jobs are the files we want to keep
        """
        return []

    def _commands_to_shell_script(self, commands, fp, shebang='#!/bin/bash'):
        f = open(fp, 'w')
        f.write(shebang)
        f.write('\n')
        f.write('\n'.join(commands))
        f.write('\n')
        f.close()


class ParallelBetaDiversitySingle(ParallelBetaDiversity):

    def _identify_files_to_remove(self, job_result_filepaths, params):
        """ The output of the individual jobs are the files we want to keep
        """
        return job_result_filepaths

    def _write_merge_map_file(self,
                              input_file_basename,
                              job_result_filepaths,
                              params,
                              output_dir,
                              merge_map_filepath):

        merge_map_f = open(merge_map_filepath, 'w')

        for metric in params['metrics'].split(','):
            fps_to_merge = [
                fp for fp in job_result_filepaths if '/%s_' %
                metric in fp]
            output_fp = join(
                output_dir, '%s_%s.txt' %
                (metric, input_file_basename))
            merge_map_f.write(
                '%s\t%s\n' %
                ('\t'.join(fps_to_merge), output_fp))

        merge_map_f.close()

    def _get_job_commands(self,
                          input_fp,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate beta diversity to split single OTU table to multiple jobs

        full_tree=True is faster: beta_diversity.py -f will make things
        go faster, but be sure you already have the correct minimal tree.
        """
        commands = []
        result_filepaths = []

        sids = load_table(input_fp).ids()

        if params['full_tree']:
            full_tree_str = '-f'
        else:
            full_tree_str = ''

        if params['tree_path']:
            tree_str = '-t %s' % params['tree_path']
        else:
            tree_str = ''

        metrics = params['metrics']

        # this is a little bit of an abuse of _merge_to_n_commands, so may
        # be worth generalizing that method - this determines the correct
        # number of samples to process in each command
        sample_id_groups = self._merge_to_n_commands(sids,
                                                     params['jobs_to_start'],
                                                     delimiter=',',
                                                     command_prefix='',
                                                     command_suffix='')

        for i, sample_id_group in enumerate(sample_id_groups):
            working_dir_i = join(working_dir, str(i))
            create_dir(working_dir_i)
            output_dir_i = join(output_dir, str(i))
            create_dir(output_dir_i)
            result_filepaths.append(output_dir_i)
            input_dir, input_fn = split(input_fp)
            input_basename, input_ext = splitext(input_fn)
            sample_id_desc = sample_id_group.replace(',', '_')
            output_fns = ['%s_%s.txt' % (metric, input_basename)
                          for metric in metrics.split(',')]
            rename_command, current_result_filepaths = self._get_rename_command(
                output_fns, working_dir_i, output_dir_i)

            result_filepaths += current_result_filepaths

            bdiv_command = '%s -i %s -o %s %s -m %s %s -r %s' %\
                (self._script_name,
                 input_fp,
                 working_dir_i,
                 tree_str,
                 params['metrics'],
                 full_tree_str,
                 sample_id_group)

            shell_script_fp = '%s/%s%d.sh' % (working_dir_i, job_prefix, i)
            shell_script_commands = [bdiv_command] + rename_command.split(';')
            self._commands_to_shell_script(shell_script_commands,
                                           shell_script_fp)
            commands.append('bash %s' % shell_script_fp)

        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)

        return commands, result_filepaths

    def _get_poller_command(self,
                            expected_files_filepath,
                            merge_map_filepath,
                            deletion_list_filepath,
                            command_prefix='/bin/bash; ',
                            command_suffix='; exit'):
        """Generate command to initiate a poller to monitior/process completed runs
        """
        result = '%s poller.py -f %s -m %s -d %s -t %d -p %s %s' % \
            (command_prefix,
             expected_files_filepath,
             merge_map_filepath,
             deletion_list_filepath,
             self._seconds_to_sleep,
             'qiime.parallel.beta_diversity.parallel_beta_diversity_process_run_results_f',
             command_suffix)

        return result, []


class ParallelBetaDiversityMultiple(ParallelBetaDiversity):

    def _get_job_commands(self,
                          input_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate beta diversity to split multiple OTU tables to multiple jobs
        """

        if params['full_tree']:
            full_tree_str = '-f'
        else:
            full_tree_str = ''

        if params['tree_path']:
            tree_str = '-t %s' % params['tree_path']
        else:
            tree_str = ''

        commands = []
        result_filepaths = []

        metrics = params['metrics']

        for input_fp in input_fps:
            input_path, input_fn = split(input_fp)
            input_basename, input_ext = splitext(input_fn)
            output_fns = \
                ['%s_%s.txt' % (metric, input_basename)
                 for metric in metrics.split(',')]
            rename_command, current_result_filepaths = self._get_rename_command(
                output_fns, working_dir, output_dir)
            result_filepaths += current_result_filepaths

            command = '%s %s -i %s -o %s %s -m %s %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 input_fp,
                 working_dir,
                 tree_str,
                 params['metrics'],
                 full_tree_str,
                 rename_command,
                 command_suffix)

            commands.append(command)

        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)

        return commands, result_filepaths


def parallel_beta_diversity_process_run_results_f(f):
    """ Handles re-assembling of a distance matrix from component vectors
    """
    # iterate over component, output fp lines
    for line in f:
        fields = line.strip().split('\t')
        dm_components = fields[:-1]
        output_fp = fields[-1]
        # assemble the current dm
        dm = assemble_distance_matrix(map(open, dm_components))
        # and write it to file
        output_f = open(output_fp, 'w')
        output_f.write(dm)
        output_f.close()

    return True


def assemble_distance_matrix(dm_components):
    """ assemble distance matrix components into a complete dm string

    """
    data = {}
    # iterate over compenents
    for c in dm_components:
        # create a blank list to store the column ids
        col_ids = []
        # iterate over lines
        for line in c:
            # split on tabs remove leading and trailing whitespace
            fields = line.strip().split()
            if fields:
                # if no column ids seen yet, these are them
                if not col_ids:
                    col_ids = fields
                # otherwise this is a data row so add it to data
                else:
                    sid = fields[0]
                    data[sid] = dict(zip(col_ids, fields[1:]))

    # grab the col/row ids as a list so it's ordered
    labels = data.keys()
    # create an empty list to build the dm
    dm = []
    # construct the dm one row at a time
    for l1 in labels:
        dm.append([float(data[l1][l2]) for l2 in labels])
    # create the dm string and return it
    dm = format_distance_matrix(labels, dm)
    return dm
