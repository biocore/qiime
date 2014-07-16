#!/usr/bin/env python
# File created on 13 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import join, split, splitext, abspath, exists, basename
from os import makedirs
from tempfile import mkdtemp
from math import ceil

import networkx as nx
from skbio.util.misc import create_dir
from biom import load_table

from qiime.parallel.util import ParallelWrapper, merge_files_from_dirs
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp
from qiime.format import format_distance_matrix


def merge_distance_matrix(output_fp, component_files):
    """Merges the distance matrix stored in multiple components on files"""
    data = {}
    # Iterate over component files
    for comp_file in component_files:
        # Create a blank list to store the columns ids
        col_ids = []
        with open(comp_file, 'U') as comp_f:
            # Iterate over lines
            for line in comp_f:
                # Split on tabs removing leading and trailing whitespace
                fields = line.strip().split()
                if fields:
                    # If no column ids seen yet, these are them
                    if not col_ids:
                        col_ids = fields
                    # Otherwise this is a data row so add it to data
                    else:
                        sid = fields[0]
                        data[sid] = dict(zip(col_ids, fields[1:]))
    # Grab the col/row ids as a lit so it's ordered
    labels = data.keys()
    # Create an empty list to build the dm
    dm = []
    # construct the dm one row at a time
    for l1 in labels:
        dm.append([data[l1][l2] for l2 in labels])
    # Store the distance matrix string
    with open(output_fp, 'w') as out_f:
        out_f.write(format_distance_matrix(labels, dm))


class ParallelBetaDiversitySingle(ParallelWrapper):
    def _get_sample_id_groups(self, biom_fp, num_groups):
        sample_ids = load_table(biom_fp).ids()
        ids_per_group = int(ceil(len(sample_ids)/float(num_groups)))
        sample_id_groups = []
        start = 0
        end = ids_per_group
        for i in range(num_groups):
            sample_id_groups.append(sample_ids[start:end])
            start = end
            end += ids_per_group
        return sample_id_groups

    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        # Create the output directory
        if not exists(output_dir):
            makedirs(output_dir)

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Create a working directory
        working_dir = mkdtemp(prefix='beta_div_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Generate the log file
        self._log_file = generate_log_fp(output_dir)

        full_tree_str = '-f' if params['full_tree'] else ''
        tree_str = ('-t %s' % abspath(params['tree_path'])
                    if params['tree_path'] else '')
        metrics = params['metrics'].split(',')

        # Determine how many samples will be computed by each command
        sample_id_groups = self._get_sample_id_groups(input_fp, jobs_to_start)

        # Construct the commands
        output_dirs = []
        node_names = []
        for i, sample_id_group in enumerate(sample_id_groups):
            node_name = "BDIV_%d" % i
            node_names.append(node_name)
            out_dir = join(working_dir, node_name)
            output_dirs.append(out_dir)

            cmd = str("beta_diversity.py -i %s -o %s %s -m %s %s -r %s"
                      % (input_fp, out_dir, tree_str, params['metrics'],
                         full_tree_str, ",".join(sample_id_group)))
            self._job_graph.add_node(node_name, job=(cmd, ),
                                     requires_deps=False)

        # Merge the results
        merge_nodes = []
        prefix = splitext(basename(input_fp))[0]
        for metric in metrics:
            merge_node = "MERGE_%s" % metric
            merge_nodes.append(merge_node)
            output_fp = join(output_dir, "%s_%s.txt" % (metric, prefix))
            format_str = "%s_*.txt" % metric
            self._job_graph.add_node(merge_node,
                                     job=(merge_files_from_dirs,
                                          output_fp, output_dirs, format_str,
                                          merge_distance_matrix),
                                     requires_deps=False)
        # Make sure that the merge nodes are executed after all the worker
        # nodes are done
        for merge_node in merge_nodes:
            for worker_node in node_names:
                self._job_graph.add_edge(worker_node, merge_node)


class ParallelBetaDiversityMultiple(ParallelWrapper):
    pass


# class ParallelBetaDiversity(ParallelWrapper):
#     _script_name = "beta_diversity.py"
#     _input_splitter = ParallelWrapper._input_existing_filepaths
#     _job_prefix = 'BDIV'

#     def _identify_files_to_remove(self, job_result_filepaths, params):
#         """ The output of the individual jobs are the files we want to keep
#         """
#         return []

#     def _commands_to_shell_script(self, commands, fp, shebang='#!/bin/bash'):
#         f = open(fp, 'w')
#         f.write(shebang)
#         f.write('\n')
#         f.write('\n'.join(commands))
#         f.write('\n')
#         f.close()


# class ParallelBetaDiversitySingle(ParallelBetaDiversity):

#     def _identify_files_to_remove(self, job_result_filepaths, params):
#         """ The output of the individual jobs are the files we want to keep
#         """
#         return job_result_filepaths

#     def _write_merge_map_file(self,
#                               input_file_basename,
#                               job_result_filepaths,
#                               params,
#                               output_dir,
#                               merge_map_filepath):

#         merge_map_f = open(merge_map_filepath, 'w')

#         for metric in params['metrics'].split(','):
#             fps_to_merge = [
#                 fp for fp in job_result_filepaths if '/%s_' %
#                 metric in fp]
#             output_fp = join(
#                 output_dir, '%s_%s.txt' %
#                 (metric, input_file_basename))
#             merge_map_f.write(
#                 '%s\t%s\n' %
#                 ('\t'.join(fps_to_merge), output_fp))

#         merge_map_f.close()

#     def _get_job_commands(self,
#                           input_fp,
#                           output_dir,
#                           params,
#                           job_prefix,
#                           working_dir,
#                           command_prefix='/bin/bash; ',
#                           command_suffix='; exit'):
#         """Generate beta diversity to split single OTU table to multiple jobs

#         full_tree=True is faster: beta_diversity.py -f will make things
#         go faster, but be sure you already have the correct minimal tree.
#         """
#         commands = []
#         result_filepaths = []

#         sids = load_table(input_fp).ids()

#         if params['full_tree']:
#             full_tree_str = '-f'
#         else:
#             full_tree_str = ''

#         if params['tree_path']:
#             tree_str = '-t %s' % params['tree_path']
#         else:
#             tree_str = ''

#         metrics = params['metrics']

#         # this is a little bit of an abuse of _merge_to_n_commands, so may
#         # be worth generalizing that method - this determines the correct
#         # number of samples to process in each command
#         sample_id_groups = self._merge_to_n_commands(sids,
#                                                      params['jobs_to_start'],
#                                                      delimiter=',',
#                                                      command_prefix='',
#                                                      command_suffix='')

#         for i, sample_id_group in enumerate(sample_id_groups):
#             working_dir_i = join(working_dir, str(i))
#             create_dir(working_dir_i)
#             output_dir_i = join(output_dir, str(i))
#             create_dir(output_dir_i)
#             result_filepaths.append(output_dir_i)
#             input_dir, input_fn = split(input_fp)
#             input_basename, input_ext = splitext(input_fn)
#             sample_id_desc = sample_id_group.replace(',', '_')
#             output_fns = ['%s_%s.txt' % (metric, input_basename)
#                           for metric in metrics.split(',')]
#             rename_command, current_result_filepaths = self._get_rename_command(
#                 output_fns, working_dir_i, output_dir_i)

#             result_filepaths += current_result_filepaths

#             bdiv_command = '%s -i %s -o %s %s -m %s %s -r %s' %\
#                 (self._script_name,
#                  input_fp,
#                  working_dir_i,
#                  tree_str,
#                  params['metrics'],
#                  full_tree_str,
#                  sample_id_group)

#             shell_script_fp = '%s/%s%d.sh' % (working_dir_i, job_prefix, i)
#             shell_script_commands = [bdiv_command] + rename_command.split(';')
#             self._commands_to_shell_script(shell_script_commands,
#                                            shell_script_fp)
#             commands.append('bash %s' % shell_script_fp)

#         commands = self._merge_to_n_commands(commands,
#                                              params['jobs_to_start'],
#                                              command_prefix=command_prefix,
#                                              command_suffix=command_suffix)

#         return commands, result_filepaths

#     def _get_poller_command(self,
#                             expected_files_filepath,
#                             merge_map_filepath,
#                             deletion_list_filepath,
#                             command_prefix='/bin/bash; ',
#                             command_suffix='; exit'):
#         """Generate command to initiate a poller to monitior/process completed runs
#         """
#         result = '%s poller.py -f %s -m %s -d %s -t %d -p %s %s' % \
#             (command_prefix,
#              expected_files_filepath,
#              merge_map_filepath,
#              deletion_list_filepath,
#              self._seconds_to_sleep,
#              'qiime.parallel.beta_diversity.parallel_beta_diversity_process_run_results_f',
#              command_suffix)

#         return result, []


# class ParallelBetaDiversityMultiple(ParallelBetaDiversity):

#     def _get_job_commands(self,
#                           input_fps,
#                           output_dir,
#                           params,
#                           job_prefix,
#                           working_dir,
#                           command_prefix='/bin/bash; ',
#                           command_suffix='; exit'):
#         """Generate beta diversity to split multiple OTU tables to multiple jobs
#         """

#         if params['full_tree']:
#             full_tree_str = '-f'
#         else:
#             full_tree_str = ''

#         if params['tree_path']:
#             tree_str = '-t %s' % params['tree_path']
#         else:
#             tree_str = ''

#         commands = []
#         result_filepaths = []

#         metrics = params['metrics']

#         for input_fp in input_fps:
#             input_path, input_fn = split(input_fp)
#             input_basename, input_ext = splitext(input_fn)
#             output_fns = \
#                 ['%s_%s.txt' % (metric, input_basename)
#                  for metric in metrics.split(',')]
#             rename_command, current_result_filepaths = self._get_rename_command(
#                 output_fns, working_dir, output_dir)
#             result_filepaths += current_result_filepaths

#             command = '%s %s -i %s -o %s %s -m %s %s %s %s' %\
#                 (command_prefix,
#                  self._script_name,
#                  input_fp,
#                  working_dir,
#                  tree_str,
#                  params['metrics'],
#                  full_tree_str,
#                  rename_command,
#                  command_suffix)

#             commands.append(command)

#         commands = self._merge_to_n_commands(commands,
#                                              params['jobs_to_start'],
#                                              command_prefix=command_prefix,
#                                              command_suffix=command_suffix)

#         return commands, result_filepaths


# def parallel_beta_diversity_process_run_results_f(f):
#     """ Handles re-assembling of a distance matrix from component vectors
#     """
#     # iterate over component, output fp lines
#     for line in f:
#         fields = line.strip().split('\t')
#         dm_components = fields[:-1]
#         output_fp = fields[-1]
#         # assemble the current dm
#         dm = assemble_distance_matrix(map(open, dm_components))
#         # and write it to file
#         output_f = open(output_fp, 'w')
#         output_f.write(dm)
#         output_f.close()

#     return True


# def assemble_distance_matrix(dm_components):
#     """ assemble distance matrix components into a complete dm string

#     """
#     print "I get called."
#     data = {}
#     # iterate over compenents
#     for c in dm_components:
#         # create a blank list to store the column ids
#         col_ids = []
#         # iterate over lines
#         for line in c:
#             # split on tabs remove leading and trailing whitespace
#             fields = line.strip().split()
#             if fields:
#                 # if no column ids seen yet, these are them
#                 if not col_ids:
#                     col_ids = fields
#                 # otherwise this is a data row so add it to data
#                 else:
#                     sid = fields[0]
#                     data[sid] = dict(zip(col_ids, fields[1:]))

#     # grab the col/row ids as a list so it's ordered
#     labels = data.keys()
#     # create an empty list to build the dm
#     dm = []
#     # construct the dm one row at a time
#     for l1 in labels:
#         dm.append([data[l1][l2] for l2 in labels])
#     # create the dm string and return it
#     dm = format_distance_matrix(labels, dm)
#     return dm
