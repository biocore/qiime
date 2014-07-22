#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import basename, join, abspath, exists, splitext
from os import makedirs
from tempfile import mkdtemp

import networkx as nx
from brokit.formatdb import build_blast_db_from_fasta_path
from skbio.parse.sequences import parse_fasta

from qiime.parallel.context import context
from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import (input_fasta_splitter, merge_files_from_dirs,
                                 concatenate_files, BufferedWriter)
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


def merge_otu_maps(output_fp, otu_maps):
    """"""
    unique_otu_map = {}
    for fp in otu_maps:
        with open(fp, 'U') as otu_map:
            for line in otu_map:
                fields = line.strip().split()
                try:
                    # current otu_id already exists, so append this set
                    # of seq_ids
                    unique_otu_map[fields[0]].extend(fields[1:])
                except KeyError:
                    # Current otu_id has not been seen yet, so create it with
                    # current set of otus
                    unique_otu_map[fields[0]] = fields[1:]

    with open(output_fp, 'w') as outf:
        for otu_id, seq_ids in unique_otu_map.items():
            outf.write('\t'.join([otu_id] + seq_ids))
            outf.write('\n')


class ParallelPickOtus(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        # Create the output directory if it does not exists
        output_dir = abspath(output_dir)
        if not exists(output_dir):
            makedirs(output_dir)

        # Generate the log file
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='otu_picker_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Perform any work that the specific OTU picker needs to do and get
        # back the list of node names that the PPOTU_X jobs should wait for,
        # and a boolean to know if the command wrapper should look for the
        # blast db or not
        dep_job_names, needs_blast = self._picker_specific_nodes(params,
                                                                 working_dir)

        # Split the input fasta file'
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)
        dep_job_names.append("SPLIT_FASTA")

        # Build the commands
        output_dirs = []
        node_names = []
        picker_specific_param_str = self._get_specific_params_str(params)
        for i in range(jobs_to_start):
            node_name = "PPOTU_%s" % i
            node_names.append(node_name)
            out_dir = join(working_dir, node_name)
            output_dirs.append(out_dir)
            cmd = ("pick_otus.py %s -o %s %s"
                   % ("-i %s", out_dir, picker_specific_param_str))
            self._job_graph.add_node(node_name,
                                     job=(command_wrapper, cmd, i,
                                          needs_blast),
                                     requires_deps=True)
            for node in dep_job_names:
                self._job_graph.add_edge(node, node_name)

        # Generate the paths to the output files
        prefix = splitext(basename(input_fp))[0]
        out_log = join(output_dir, "%s_otus.log" % prefix)
        out_otus = join(output_dir, "%s_otus.txt" % prefix)
        out_failures = join(output_dir, "%s_failures.txt" % prefix)
        self._job_graph.add_node("MERGE_LOGS",
                                 job=(merge_files_from_dirs, out_log,
                                      output_dirs, "*_otus.log",
                                      concatenate_files),
                                 requires_deps=False)
        self._job_graph.add_node("MERGE_FAILURES",
                                 job=(merge_files_from_dirs, out_failures,
                                      output_dirs, "*_failures.txt",
                                      concatenate_files),
                                 requires_deps=False)
        self._job_graph.add_node("MERGE_OTU_MAPS",
                                 job=(merge_files_from_dirs, out_otus,
                                      output_dirs, "*_otus.txt",
                                      merge_otu_maps),
                                 requires_deps=False)
        # Make sure that the merge nodes are executed after the worker nodes
        for node in node_names:
            self._job_graph.add_edge(node, "MERGE_LOGS")
            self._job_graph.add_edge(node, "MERGE_FAILURES")
            self._job_graph.add_edge(node, "MERGE_OTU_MAPS")


class ParallelPickOtusUclustRef(ParallelPickOtus):
    def _picker_specific_nodes(self, params, working_dir):
        return [], False

    def _get_specific_params_str(self, params):
        enable_rev_strand_match_str = (
            '-z' if params['enable_rev_strand_match'] else '')
        optimal_uclust_str = '-A' if params['optimal_uclust'] else ''
        exact_uclust_str = '-E' if params['exact_uclust'] else ''
        stable_sort_str = (
            '' if params['stable_sort'] else '--suppress_uclust_stable_sort')
        save_uc_files_str = '' if params['save_uc_files'] else '-d'

        return (
            "-r %s -m uclust_ref --suppress_new_clusters -s %s %s %s %s "
            "--max_accepts %s --max_rejects %s --stepwords %d --w %d %s %s"
            % (params['refseqs_fp'], params['similarity'],
               enable_rev_strand_match_str, optimal_uclust_str,
               exact_uclust_str, params['max_accepts'], params['max_rejects'],
               params['stepwords'], params['word_length'], stable_sort_str,
               save_uc_files_str))


class ParallelPickOtusUsearch61Ref(ParallelPickOtus):
    def _picker_specific_nodes(self, params, working_dir):
        return [], False

    def _get_specific_params_str(self, params):
        # Generate the parameters to pass to pick_otus.py. This must exclude
        # parameters that get passed only to the parallel version
        # (e.g jobs_to_start) and values that get overwritten (e.g.,
        # input_fasta_fp)
        param_fields = []
        ignored_params = {"input_fasta_fp", "output_dir", "jobs_to_start",
                          "retain_temp_files", "suppress_submit_jobs",
                          "poll_directly", "cluster_jobs_fp",
                          "suppress_polling", "job_prefix", "seconds_to_sleep"}
        for name, value in params.items():
            if name in ignored_params or value is False:
                pass
            elif value is True:
                param_fields.append('--%s' % name)
            else:
                param_fields.append('--%s %s' % (name, value))
        return ("-m usearch61_ref --suppress_new_clusters %s"
                % ' '.join(param_fields))


class ParallelPickOtusBlast(ParallelPickOtus):
    def _picker_specific_nodes(self, params, working_dir):
        node = "BUILD_BLAST_DB"
        if not params['blast_db']:
            # Build the blast database from the refseqs_fp -- all procs
            # will then access one db rather than create one per proc
            self._job_graph.add_node("BUILD_BLAST_DB",
                                     job=(build_blast_db_from_fasta_path,
                                          params['refseqs_fp'], False,
                                          working_dir, False),
                                     requires_deps=False)
        return [node], True

    def _get_specific_params_str(self, params):
        return (
            "-m blast -e %s -s %s --min_aligned_percent %s %s"
            % (params['max_e_value'], params['similarity'],
               params['min_aligned_percent'], "-b %s"))


class ParallelPickOtusTrie(ParallelWrapper):
    """Picking Otus using a trie the parallel way

    We parallelize the Trie OTU picker using this scheme:
    1. use the exact prefix filter with a short wordlength (say 5)
       to sort all reads into buckets according to their first 5 nucs.
    2. Run trie otupicker an each bucket separately, distribute over cluster
    3. Combine mappings of 2, Since each bucket is independent fro the rest,
       a simple cat with incrementing OTU ids should do it.

    Since the parallel trie scheme is different from the other pick otu
    algorithms, we are not extending the ParallelPickOtus class
    """
    def _construct_job_graph(self, input_fp, output_dir, params):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        # Create the output directory if it does not exists
        output_dir = abspath(output_dir)
        if not exists(output_dir):
            makedirs(output_dir)

        # Generate the log file
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='otu_picker_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Split the input fasta file. We cannot add the job to the job graph
        # because the number of jobs will be defined by the number of files
        # that the input file is split
        prefix_length = params['prefix_length'] or 1
        fasta_fps = trie_input_splitter(input_fp, working_dir, prefix_length)

        # Build the commands
        output_dirs = []
        node_names = []
        for i, fasta_fp in enumerate(fasta_fps):
            node_name = "PPOTU_%s" % i
            node_names.append(node_name)
            out_dir = join(working_dir, node_name)
            output_dirs.append(out_dir)
            cmd = ("pick_otus.py -i %s -o %s -m trie"
                   % (fasta_fp, out_dir))
            self._job_graph.add_node(node_name, job=(cmd,),
                                     requires_deps=False)

        # Generate the paths to the output files
        prefix = splitext(basename(input_fp))[0]
        out_log = join(output_dir, "%s_otus.log" % prefix)
        out_otus = join(output_dir, "%s_otus.txt" % prefix)
        self._job_graph.add_node("MERGE_LOGS",
                                 job=(merge_files_from_dirs, out_log,
                                      output_dirs, "*_otus.log",
                                      concatenate_files),
                                 requires_deps=False)
        self._job_graph.add_node("MERGE_OTU_MAPS",
                                 job=(merge_files_from_dirs, out_otus,
                                      output_dirs, "*_otus.txt",
                                      merge_otu_maps_trie),
                                 requires_deps=False)
        # Make sure that the merge nodes are executed after the worker nodes
        for node in node_names:
            self._job_graph.add_edge(node, "MERGE_LOGS")
            self._job_graph.add_edge(node, "MERGE_OTU_MAPS")


def trie_input_splitter(fasta_fp, output_dir, prefix_length):
    """Split input sequences into sets with identical prefix

    Parameters
    ----------
    fasta_fp : str
        Path to the input fasta file
    output_dir : str
        Path to the output directory
    prefix_length : int
        The length of the prefix
    """
    out_files = []
    buffered_handles = {}
    with open(fasta_fp, 'U') as f:
        for seq_id, seq in parse_fasta(f):
            if len(seq) < prefix_length:
                raise ValueError("Prefix length must be equal or shorter than "
                                 "than sequence.\n Found seq %s with length "
                                 "%d" % (seq_id, len(seq)))
            prefix = seq[:prefix_length]

            if prefix not in buffered_handles:
                # Never seen this prefix before
                out_fp = join(output_dir, prefix)
                buffered_handles[prefix] = BufferedWriter(out_fp)
                out_files.append(out_fp)

            buffered_handles[prefix].write('>%s\n%s\n' % (seq_id, seq))

    # make sure all buffers are closed and flushed
    for buf_fh in buffered_handles.itervalues():
        buf_fh.close()

    return out_files


def merge_otu_maps_trie(output_fp, otu_maps):
    """"""
    otu_id = 0
    with open(output_fp, 'w') as outf:
        for fp in otu_maps:
            with open(fp, 'U') as otu_map:
                for line in otu_map:
                    fields = line.strip().split()
                    outf.write('\t'.join(["%d" % otu_id] + fields[1:]))
                    outf.write('\n')
                    otu_id += 1
