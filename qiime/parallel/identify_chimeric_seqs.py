#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import abspath, basename, join, exists
from os import makedirs
from shutil import copyfile
from tempfile import mkdtemp

import networkx as nx
from brokit.formatdb import build_blast_db_from_fasta_path
from skbio.parse.sequences import parse_fasta

from qiime.identify_chimeric_seqs import make_cidx_file
from qiime.util import write_degapped_fasta_to_file
from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import (input_fasta_splitter, concatenate_files,
                                 command_wrapper, fasta_splitter_handler,
                                 blast_db_builder_handler)
from qiime.parallel.context import context, system_call
from qiime.workflow.util import generate_log_fp, WorkflowLogger


class ParallelChimericSeqIdentifier(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        """Builds the job workflow graph"""
        # Do the parameter parsing
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        method = params['chimera_detection_method']

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Create the output directory
        if not exists(output_dir):
            makedirs(output_dir)

        # Generate the log
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='ICS_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        dep_job_names, keys, funcs = self._csi_specific_nodes(params,
                                                              working_dir)

        # Split the input_fp in multiple files
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)
        dep_job_names.append("SPLIT_FASTA")
        keys.insert(0, "SPLIT_FASTA")
        funcs["SPLIT_FASTA"] = fasta_splitter_handler

        # Create the commands
        temp_out_fps = []
        ics_nodes = []
        method_parms = self._get_specific_params_str(params)
        for i in range(jobs_to_start):
            temp_out_fp = join(working_dir, "chimeric_seqs_%d.txt" % i)
            temp_out_fps.append(temp_out_fp)
            node_name = "ICS_%d" % i
            ics_nodes.append(node_name)
            cmd = ("identify_chimeric_seqs.py -o %s %s %s"
                   % (temp_out_fp, "-i %s", method_parms))
            self._job_graph.add_node(node_name,
                                     job=(command_wrapper, cmd, i, keys,
                                          funcs),
                                     requires_deps=True)
            # Adding the dependency edges to the graph
            for dep_node_name in dep_job_names:
                self._job_graph.add_edge(dep_node_name, node_name)

        # Adding temp_out_fps to the clean_up variable
        self._filepaths_to_remove.extend(temp_out_fps)
        # Merge the results by concatenating the output files
        output_fp = join(output_dir, "chimeric_seqs.txt")
        self._job_graph.add_node("CONCAT", job=(concatenate_files,
                                                output_fp, temp_out_fps),
                                 requires_deps=False)
        # Make sure that the "CONCAT" node is the latest executed node
        for node in ics_nodes:
            self._job_graph.add_edge(node, "CONCAT")


class ParallelChimericSeqIdentifierBlast(ParallelChimericSeqIdentifier):
    def _csi_specific_nodes(self, params, working_dir):
        node = "BUILD_BLAST_DB"
        # Build the blast database
        self._job_graph.add_node(node,
                                 job=(build_blast_db_from_fasta_path,
                                      params['reference_seqs_fp'], False,
                                      working_dir, False),
                                 requires_deps=False)
        return [node], [node], {node: blast_db_builder_handler}

    def _get_specific_params_str(self, params):
        return ("-t {0} -m blast_fragments -n {1} -d {2} -e {3} -b %s".format(
            params['id_to_taxonomy_fp'], params['num_fragments'],
            params['taxonomy_depth'], params['max_e_value']))


class ParallelChimericSeqIdentifierChimSlayer(ParallelChimericSeqIdentifier):
    def _csi_specific_nodes(self, params, working_dir):
        # In the ChimeraSlayer method, we have to perform numerous
        # steps prior to the actual execution of the jobs.
        # We first need to copy the reference files to the working
        # directory because ChimeraSlayer creates an index file of the
        # reference and it will crash without write permission in the
        # directory in which the reference sequences file is
        aln_ref_seqs_fp = abspath(params['aligned_reference_seqs_fp'])
        dest_aln_ref_seqs_fp = join(working_dir, basename(aln_ref_seqs_fp))
        copyfile(aln_ref_seqs_fp, dest_aln_ref_seqs_fp)
        aln_ref_seqs_fp = dest_aln_ref_seqs_fp
        # Adding the copy to the clean up variable
        self._filepaths_to_remove.append(aln_ref_seqs_fp)
        ref_seqs_fp = params['reference_seqs_fp']

        if ref_seqs_fp:
            ref_seqs_fp = abspath(ref_seqs_fp)
            # Copy the unaligned reference database
            dest_ref_seqs_fp = join(working_dir, basename(ref_seqs_fp))
            copyfile(ref_seqs_fp, dest_ref_seqs_fp)
            ref_seqs_fp = dest_ref_seqs_fp
        else:
            # The unaligned reference database does not exists, create it
            with open(aln_ref_seqs_fp, 'U') as f:
                ref_seqs_fp = write_degapped_fasta_to_file(
                    parse_fasta(f), tmp_dir=working_dir)
        # Add the copied/new ref_seqs_fp to the clean up variable
        self._filepaths_to_remove.append(ref_seqs_fp)

        # Make the index file. ChimeraSlayer first checks to see if the
        # index file exists. If not, it tries to create it. This can lead
        # to race conditions if several parallel jobs try to create it at
        # the same time
        self._job_graph.add_node("MAKE_CIDX", job=(make_cidx_file,
                                                   aln_ref_seqs_fp),
                                 requires_deps=False)
        self._filepaths_to_remove.append("%s.cidx" % aln_ref_seqs_fp)
        params['aligned_reference_seqs_fp'] = aln_ref_seqs_fp
        params['reference_seqs_fp'] = ref_seqs_fp

        return ["MAKE_CIDX"], [], {}

    def _get_specific_params_str(self, params):
        min_div_ratio_str = ("--min_div_ratio %s" % params['min_div_ratio']
                             if params['min_div_ratio'] else "")
        return (" -m ChimeraSlayer -a %s -r %s %s"
                % (params['aligned_reference_seqs_fp'],
                   params['reference_seqs_fp'],
                   min_div_ratio_str))
