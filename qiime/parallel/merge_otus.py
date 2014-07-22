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

# from os.path import basename, join
# from time import time
# from cogent.core.tree import TreeNode
# from os import system
# import os


# class JobError(Exception):
#     pass

# INTERNAL_COUNT = 0


# def mergetree(left, right, working_dir):
#     """Reconstruct a tree from merge order"""
#     # decorate and infer filenames for tips
#     global INTERNAL_COUNT

#     if not isinstance(left, TreeNode):
#         filepath = str(left[0])
#         name = basename(filepath.split('.')[0])
#         left = TreeNode(Name=name)
#         left.FilePath = filepath
#         left.Processed = False
#         left.PollPath = None  # doesn't make sense for tips
#         left.FullCommand = None
#         left.EndTime = None
#         left.StartTime = None
#         left.TotalTime = None

#     if not isinstance(right, TreeNode):
#         filepath = str(right[0])
#         name = basename(filepath.split('.')[0])
#         right = TreeNode(Name=name)
#         right.FilePath = filepath
#         right.Processed = False
#         right.PollPath = None  # doesn't make sense for tips
#         right.FullCommand = None
#         right.EndTime = None
#         right.StartTime = None
#         right.TotalTime = None

#     # internal node
#     name = str(INTERNAL_COUNT)
#     filepath = join(working_dir, name) + '.biom'
#     merged = TreeNode(Name=name, Children=[left, right])
#     merged.FilePath = filepath
#     merged.Processed = False
#     merged.PollPath = filepath + '.poll'
#     merged.FullCommand = None
#     merged.EndTime = None
#     merged.StartTime = None
#     merged.TotalTime = None

#     INTERNAL_COUNT += 1
#     return merged


# def reset_internal_count():
#     global INTERNAL_COUNT
#     INTERNAL_COUNT = 0


# def mergeorder(items, working_dir):
#     """Code taken from http://en.literateprograms.org/Merge_sort_(Python)"""
#     if len(items) < 2:
#         return items
#     middle = len(items) / 2
#     left = mergeorder(items[:middle], working_dir)
#     right = mergeorder(items[middle:], working_dir)
#     return mergetree(left, right, working_dir)


# def initial_nodes_to_merge(tree):
#     """Determine what nodes are safe to process first

#     The first nodes to process are those internal nodes that have tips as
#     children
#     """
#     to_process = set([])
#     for n in tree.tips():
#         sibs_are_tips = [s.istip() for s in n.siblings()]
#         if all(sibs_are_tips):
#             to_process.add(n.Parent)
#     return to_process


# def initial_has_dependencies(tree, to_process):
#     """All nodes that aren't processed up front stll have dependencies

#     to_process : set of nodes that are in the first round of processing
#     """
#     has_dependencies = []
#     for n in tree.nontips(include_self=True):
#         if n not in to_process:
#             has_dependencies.append(n)
#     return has_dependencies


# def job_complete(node, verbose=False):
#     """Check if the job is complete"""
#     if node.PollPath is None or node.istip():
#         raise JobError("Attempting to merge tip: %s" % node.Name)

#     if node.Processed:
#         raise JobError("Already processed node: %s" % node.Name)

#     if os.path.exists(node.PollPath):
#         node.EndTime = time()
#         node.TotalTime = node.EndTime - node.StartTime

#         node.ExitStatus = open(node.PollPath).read().strip()
#         if node.ExitStatus != '0':
#             raise JobError("Node %s did not complete correctly!" % node.Name)

#         if verbose:
#             print "finishing %s, %f seconds" % (node.Name, node.TotalTime)

#         node.Processed = True
#         return True

#     else:
#         return False


# def torque_job(cmd, pollpath, name, queue):
#     """Wrap a cmd for job submission"""
#     qsub_call = "qsub -k oe -N %s -q %s" % ("MOTU", queue)
#     to_submit = 'echo "%s; echo $? > %s" | %s' % (cmd, pollpath, qsub_call)

#     return to_submit


# def local_job(cmd, pollpath, name, queue):
#     """make a local job"""
#     to_submit = '%s; echo $? > %s' % (cmd, pollpath)

#     return to_submit


# def start_job(node, merge_otus_fp, queue, wrap_call=torque_job, submit=True):
#     """Starts a process"""
#     strfmt = {'MergeOTUs': merge_otus_fp,
#               'Output': node.FilePath,
#               'BIOM_A': node.Children[0].FilePath,
#               'BIOM_B': node.Children[1].FilePath}

#     cmd = "%(MergeOTUs)s -i %(BIOM_A)s,%(BIOM_B)s -o %(Output)s"
#     wrapped = wrap_call(cmd % strfmt, node.PollPath, node.Name, queue)

#     if submit:
#         system(wrapped)

#     node.FullCommand = wrapped
#     node.StartTime = time()
