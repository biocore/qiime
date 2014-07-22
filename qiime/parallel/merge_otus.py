__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from os.path import abspath, exists, splitext, basename, join
from os import makedirs, rename
from tempfile import mkdtemp

import networkx as nx

from qiime.parallel.util import ParallelWrapper
from qiime.workflow.util import generate_log_fp


class ParallelMergeOtus(ParallelWrapper):
    def _construct_job_graph(self, input_fps, output_dir):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        # Create the output directory
        output_dir = abspath(output_dir)
        if not exists(output_dir):
            makedirs(output_dir)

        # Generate the log file
        self._log_file = generate_log_fp(output_dir)

        # Create the working directory
        working_dir = mkdtemp(prefix='MOTU_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # Initiate the internal node count
        self._internal_count = 0

        # Build the tree
        node = self._mergeorder(input_fps, working_dir)

        # We need to rename the last file to the actual output file
        old_fp = self._job_graph.node[node]['output_fp']
        output_fp = join(output_dir, "merged.biom")
        self._job_graph.add_node("RENAME", job=(rename, old_fp, output_fp),
                                 requires_deps=False)

        # Make sure that the renaming happens after the final OTU table
        # have been generated
        self._job_graph.add_edge(node, "RENAME")

    def _mergetree(self, left, right, working_dir):
        """Reconstruct a tree from merge order"""
        # internal node
        nodename = str(self._internal_count)
        self._internal_count += 1

        left_fp = self._job_graph.node[left]['output_fp']
        right_fp = self._job_graph.node[right]['output_fp']
        output_fp = join(working_dir, nodename)

        cmd = "merge_otu_tables.py -i %s,%s -o %s" % (left_fp, right_fp,
                                                      output_fp)
        self._job_graph.add_node(nodename, job=(cmd,), requires_deps=False,
                                 output_fp=output_fp)
        self._job_graph.add_edge(left, nodename)
        self._job_graph.add_edge(right, nodename)
        return nodename

    def _mergeorder(self, items, working_dir):
        """Code adapted from http://en.literateprograms.org/Merge_sort_(Python)
        """
        if len(items) == 1:
            output_fp = abspath(items[0])
            nodename = splitext(basename(output_fp))[0]
            f = lambda: None
            self._job_graph.add_node(nodename, job=(f,), requires_deps=False,
                                     output_fp=output_fp)
            return nodename
        middle = len(items) / 2
        left = self._mergeorder(items[:middle], working_dir)
        right = self._mergeorder(items[middle:], working_dir)
        return self._mergetree(left, right, working_dir)
