#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import abspath, exists, basename, splitext, join
from os import makedirs
from tempfile import mkdtemp

from brokit.formatdb import build_blast_db_from_fasta_path

from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import (input_fasta_splitter, merge_files_from_dirs,
                                 concatenate_files)
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp, WorkflowLogger


def command_wrapper(cmd, idx, needs_blast, dep_results=None):
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

    if needs_blast:
        if "BUILD_BLAST_DB" not in dep_results:
            raise ValueError("Wrong job graph workflow. Node 'BUILD_BLAST_DB' "
                             "not listed as dependency of current node")
        blast_db, db_files_to_remove = dep_results["BUILD_BLAST_DB"]
        cmd = cmd % (fasta_fps[idx], blast_db)
    else:
        cmd = cmd % fasta_fps[idx]
    return system_call(cmd)


class ParallelTaxonomyAssigner(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the output directory if it does not exists
        output_dir = abspath(output_dir)
        if not exists(output_dir):
            makedirs(output_dir)

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Generate the log file
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='tax_assigner_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Perform any work that the specific tax_assigner needs to do
        # and get back the list of node names that the PTA_X jobs should wait
        # for, and if the command wrapper should look for the blast db or not
        dep_job_names, needs_blast = self._tax_specific_nodes(params,
                                                              working_dir)

        # Split the input fasta file
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)
        dep_job_names.append("SPLIT_FASTA")

        # Build the commands
        output_dirs = []
        node_names = []
        tax_specific_param_str = self._get_specific_params_str(params)
        for i in range(jobs_to_start):
            node_name = "PTA_%d" % i
            out_dir = join(working_dir, node_name)
            cmd = ("assign_taxonomy.py -o {0} -i %s {1}".format(
                out_dir, tax_specific_param_str))
            output_dirs.append(out_dir)
            node_names.append(node_name)
            self._job_graph.add_node(node_name, job=(command_wrapper, cmd, i,
                                                     needs_blast),
                                     requires_deps=True)
            # Add dependencies
            for dep_node_name in dep_job_names:
                self._job_graph.add_edge(dep_node_name, node_name)

        # Generate the paths to the output files
        prefix = splitext(basename(input_fp))[0]
        out_tax_fp = join(output_dir, "%s_tax_assignments.txt" % prefix)
        log_fp = join(output_dir, "%s_tax_assignments.log" % prefix)
        # Merge the results by concatenating the output files
        self._job_graph.add_node("CONCAT_TAX_ASSIGN",
                                 job=(merge_files_from_dirs, out_tax_fp,
                                      output_dirs, "*_tax_assignments.txt",
                                      concatenate_files),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_LOG",
                                 job=(merge_files_from_dirs, log_fp,
                                      output_dirs, "*_tax_assignments.log",
                                      concatenate_files),
                                 requires_deps=False)
        # Make sure that the concatenate jobs are executed after the worker
        # are done
        for node in node_names:
            self._job_graph.add_edge(node, "CONCAT_TAX_ASSIGN")
            self._job_graph.add_edge(node, "CONCAT_LOG")


class ParallelBlastTaxonomyAssigner(ParallelTaxonomyAssigner):
    def _tax_specific_nodes(self, params, working_dir):
        self._job_graph.add_node("BUILD_BLAST_DB",
                                 job=(build_blast_db_from_fasta_path,
                                      params['reference_seqs_fp'], False,
                                      working_dir, False),
                                 requires_deps=False)
        return ["BUILD_BLAST_DB"], True

    def _get_specific_params_str(self, params):
        return "-m blast -e {0} -b %s -t {1}".format(
            params['e_value'], params['id_to_taxonomy_fp'])


class ParallelRdpTaxonomyAssigner(ParallelTaxonomyAssigner):
    def _tax_specific_nodes(self, params, working_dir):
        return [], False

    def _get_specific_params_str(self, params):
        return ("-m rdp -c %1.2f --rdp_max_memory %d"
                % (params['confidence'], params['rdp_max_memory']))


class ParallelUclustConsensusTaxonomyAssigner(ParallelTaxonomyAssigner):
    def _tax_specific_nodes(self, params, working_dir):
        return [], False

    def _get_specific_params_str(self, params):
        return ("-m uclust --uclust_min_consensus_fraction %f "
                "--uclust_similarity %f --uclust_max_accepts %d -t %s -r %s"
                % (params['uclust_min_consensus_fraction'],
                   params['uclust_similarity'],
                   params['uclust_max_accepts'],
                   params['id_to_taxonomy_fp'],
                   params['reference_seqs_fp']))
