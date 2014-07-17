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

from qiime.parallel.context import context
from qiime.parallel.util import (ParallelWrapper, input_fasta_splitter,
                                 merge_files_from_dirs, concatenate_files)
from qiime.workflow.util import generate_log_fp


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
        self._log_file = generate_log_fp(output_dir)

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='otu_picker_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Perform any work that the specific OTU picker needs to do and get
        # back the list of node names that the PPOTU_X jobs should wait for/
        # dep_job_names = self._picker_specific_nodes(params, working_dir)

        # Split the input fasta file
        fasta_fps = input_fasta_splitter(input_fp, working_dir, jobs_to_start)

        # Build the commands
        output_dirs = []
        node_names = []
        picker_specific_param_str = self._get_specific_params_str(params)
        for i, fasta_fp in enumerate(fasta_fps):
            node_name = "PPOTU_%s" % i
            node_names.append(node_name)
            out_dir = join(working_dir, node_name)
            output_dirs.append(out_dir)
            cmd = ("pick_otus.py -i %s -o %s %s"
                   % (fasta_fp, out_dir, picker_specific_param_str))
            self._job_graph.add_node(node_name, job=(cmd,),
                                     requires_deps=False)

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
    pass


class ParallelPickOtusBlast(ParallelPickOtus):
    pass


class ParallelPickOtusTrie(ParallelPickOtus):

    """Picking Otus using a trie the parallel way

    We parallelize the Trie OTU picker using this scheme:
    1. use the exact prefix filter with a short wordlength (say 5)
       to sort all reads into buckets according to their first 5 nucs.
    2. Run trie otupicker an each bucket separately, distribute over cluster
    3. Combine mappings of 2, Since each bucket is independent fro the rest,
       a simple cat with incrementing OTU ids should do it.
    """
    pass


def greedy_partition(counts, n):
    """Distribute k counts evenly across n buckets,

    counts: dict of key, counts pairs
    n: number of buckets that the counts should be distributed over
    """

    buckets = [[] for i in range(n)]
    fill_levels = [0 for i in range(n)]

    for key in sorted(counts, reverse=True,
                      key=lambda c: counts[c]):
        smallest = fill_levels.index(min(fill_levels))
        buckets[smallest].append(key)
        fill_levels[smallest] += counts[key]

    return buckets, fill_levels
