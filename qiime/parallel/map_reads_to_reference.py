#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import abspath, exists, join
from os import makedirs
from tempfile import mkdtemp

import networkx as nx

from qiime.workflow.util import generate_log_fp
from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import (input_fasta_splitter, merge_files_from_dirs,
                                 concatenate_files)
from qiime.parallel.context import context
from qiime.parallel.pick_otus import merge_otu_maps


def command_wrapper(cmd, idx, dep_results=None):
    """Wraps the command to be executed so it can use the results produced by
    the jobs in which the command depends on

    Parameters
    ----------
    cmd : str
        Command to execute
    idx : int
        The fasta fp index that this job has to execute
    dep_results : dict of {node_name: tuple}
        The results in which cmd depends on
    """
    from qiime.parallel.context import system_call
    if "SPLIT_FASTA" not in dep_results:
        raise ValueError("Wrong job graph workflow. Node 'SPLIT_FASTA' "
                         "not listed as dependency of current node")
    fasta_fps = dep_results["SPLIT_FASTA"]

    cmd = cmd % fasta_fps[idx]

    return system_call(cmd)


def generate_biom_table(biom_fp, observation_map_fp, observation_metadata_fp):
    # Importing here so the become available on the workers
    from qiime.make_otu_table import make_otu_table
    from qiime.parse import parse_observation_metadata
    from qiime.util import write_biom_table

    # Check if we actually have observation metadata or not
    if observation_metadata_fp is not None:
        with open(abspath(observation_metadata_fp), 'U') as f:
            observation_metadata = parse_observation_metadata(f)
    else:
        observation_metadata = None

    # Create the table
    with open(observation_map_fp, 'U') as f:
        biom_table = make_otu_table(f, observation_metadata)

    # Store the table
    write_biom_table(biom_table, biom_fp)


class ParallelDatabaseMapper(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        # Create the output directory if it does not exists
        output_dir = abspath(output_dir)
        if not exists(output_dir):
            makedirs(output_dir)

        # Generate the log file
        self._log_file = generate_log_fp(output_dir)

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='database_mapper_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        refseqs_fp = abspath(params['refseqs_fp'])

        # Split the input filepath
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)

        out_dirs = []
        nodes = []
        mapper_specific_param_str = self._get_specific_params_str(params)
        for i in range(jobs_to_start):
            node = "PDM_%d" % i
            out_dir = join(working_dir, node)
            nodes.append(node)
            out_dirs.append(out_dir)
            cmd = ("map_reads_to_reference.py -i %s -r {0} -o {1} {2}".format(
                refseqs_fp, out_dir, mapper_specific_param_str))
            self._job_graph.add_node(node, job=(command_wrapper, cmd, i),
                                     requires_deps=True)
            # Make sure that any of the workers are executed until the
            # SPLIT_FASTA node is executed
            self._job_graph.add_edge("SPLIT_FASTA", node)

        # Generate the observation map
        out_obs_map = join(output_dir, 'observation_map.txt')
        self._job_graph.add_node("MERGE_OBS_MAPS",
                                 job=(merge_files_from_dirs, out_obs_map,
                                      out_dirs, "observation_map.txt",
                                      merge_otu_maps),
                                 requires_deps=False)
        # Depending on the method used, the merge of other files is done
        # different
        merge_nodes = self._get_specific_merge_nodes(out_dirs, output_dir)
        # Don't forget to add the "MERGE_OBS_MAPS" node!!
        merge_nodes.append("MERGE_OBS_MAPS")

        # Make sure that the merge commands are executed after the workers
        for node in nodes:
            for m_node in merge_nodes:
                self._job_graph.add_edge(node, m_node)

        # Create the final biom table
        biom_fp = join(output_dir, 'observation_table.biom')
        self._job_graph.add_node("MAKE_TABLE",
                                 job=(generate_biom_table, biom_fp,
                                      out_obs_map,
                                      params['observation_metadata_fp']),
                                 requires_deps=False)

        # Make sure that the creation of the OTU table happens after the
        # observation maps have been merged
        self._job_graph.add_edge("MERGE_OBS_MAPS", "MAKE_TABLE")


class ParallelDatabaseMapperUsearch(ParallelDatabaseMapper):
    def _get_specific_params_str(self, params):
        return (
            "-m usearch --min_percent_id %s --max_accepts %d --max_rejects %d "
            "--queryalnfract %f --targetalnfract %f --evalue %e"
            % (params['min_percent_id'], params['max_accepts'],
               params['max_rejects'], params['queryalnfract'],
               params['targetalnfract'], params['evalue']))

    def _get_specific_merge_nodes(self, work_dirs, output_dir):
        out_uc_fp = join(output_dir, 'out.uc')
        self._job_graph.add_node("MERGE_UCS",
                                 job=(merge_files_from_dirs, out_uc_fp,
                                      work_dirs, "out.uc", concatenate_files),
                                 requires_deps=False)
        out_bl6_fp = join(output_dir, 'out.bl6')
        self._job_graph.add_node("MERGE_BL6",
                                 job=(merge_files_from_dirs, out_bl6_fp,
                                      work_dirs, "out.bl6", concatenate_files),
                                 requires_deps=False)
        return ["MERGE_UCS", "MERGE_BL6"]


class ParallelDatabaseMapperBlat(ParallelDatabaseMapper):
    def _get_specific_params_str(self, params):
        return ("-m blat --min_percent_id %s --evalue %e"
                % (params['min_percent_id'], params['evalue']))

    def _get_specific_merge_nodes(self, work_dirs, output_dir):
        out_log_fp = join(output_dir, 'observation_table.log')
        self._job_graph.add_node("MERGE_LOGS",
                                 job=(merge_files_from_dirs, out_log_fp,
                                      work_dirs, "*.log", concatenate_files),
                                 requires_deps=False)
        out_bl9_fp = join(output_dir, 'out.bl9')
        self._job_graph.add_node("MERGE_BL9",
                                 job=(merge_files_from_dirs, out_bl9_fp,
                                      work_dirs, "*.bl9", concatenate_files),
                                 requires_deps=False)
        return ["MERGE_LOGS", "MERGE_BL9"]


class ParallelDatabaseMapperBwaShort(ParallelDatabaseMapper):
    def _get_specific_params_str(self, params):
        max_diff_str = ("--max_diff %s" % params['max_diff']
                        if params['max_diff'] is not None else "")
        return ("-m bwa-short %s " % max_diff_str)

    def _get_specific_merge_nodes(self, work_dirs, output_dir):
        out_log_fp = join(output_dir, 'observation_table.log')
        self._job_graph.add_node("MERGE_LOGS",
                                 job=(merge_files_from_dirs, out_log_fp,
                                      work_dirs, "*.log", concatenate_files),
                                 requires_deps=False)
        out_sam_fp = join(output_dir, 'bwa_raw_out.sam')
        self._job_graph.add_node("MERGE_SAM",
                                 job=(merge_files_from_dirs, out_sam_fp,
                                      work_dirs, "*.sam", concatenate_files),
                                 requires_deps=False)
        return ["MERGE_LOGS", "MERGE_SAM"]
