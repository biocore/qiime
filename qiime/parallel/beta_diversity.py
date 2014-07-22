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

from os.path import join, splitext, abspath, exists, basename
from os import makedirs
from tempfile import mkdtemp
from math import ceil
from shutil import move

import networkx as nx
from biom import load_table

from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import merge_files_from_dirs
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp, WorkflowLogger
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
        """Groups the sample ids on biom_fp in num_groups groups"""
        # Get the sample ids
        sample_ids = load_table(biom_fp).ids()
        # Compute the number of ids that has to be in each group
        ids_per_group = int(ceil(len(sample_ids)/float(num_groups)))
        # Group sample ids
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
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # Parse parameters
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
    def _construct_job_graph(self, input_fps, output_dir, params):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        output_dir = abspath(output_dir)
        # Create the output directory
        if not exists(output_dir):
            makedirs(output_dir)

        # create a working directory
        working_dir = mkdtemp(prefix='beta_div_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Generate the log file
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # Parse parameters
        full_tree_str = '-f' if params['full_tree'] else ''
        tree_str = ('-t %s' % abspath(params['tree_path'])
                    if params['tree_path'] else '')
        metrics = params['metrics'].split(',')

        # Generate the comands
        for i, input_fp in enumerate(input_fps):
            prefix = splitext(basename(input_fp))[0]
            node_name = "BDIV_%d" % i
            out_dir = join(working_dir, node_name)
            cmd = ("beta_diversity.py -i %s -o %s %s -m %s %s"
                   % (input_fp, out_dir, tree_str, params['metrics'],
                      full_tree_str))
            self._job_graph.add_node(node_name, job=(cmd,),
                                     requires_deps=False)
            # Add nodes to move the output files to the correct location
            for metric in metrics:
                fn = "%s_%s.txt" % (metric, prefix)
                cmd_fp = join(out_dir, fn)
                output_fp = join(output_dir, fn)
                move_node = "MOVE_%s_%s" % (node_name, metric)
                self._job_graph.add_node(move_node,
                                         job=(move, cmd_fp, output_fp),
                                         requires_deps=False)
                # Make sure that the move nodes are executed after the job is
                # done
                self._job_graph.add_edge(node_name, move_node)
