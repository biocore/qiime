#!/usr/bin/env python

from os.path import basename, join
from time import time
from cogent.core.tree import TreeNode
from os import system
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


class JobError(Exception):
    pass

INTERNAL_COUNT = 0


def mergetree(left, right, working_dir):
    """Reconstruct a tree from merge order"""
    # decorate and infer filenames for tips
    global INTERNAL_COUNT

    if not isinstance(left, TreeNode):
        filepath = str(left[0])
        name = basename(filepath.split('.')[0])
        left = TreeNode(Name=name)
        left.FilePath = filepath
        left.Processed = False
        left.PollPath = None  # doesn't make sense for tips
        left.FullCommand = None
        left.EndTime = None
        left.StartTime = None
        left.TotalTime = None

    if not isinstance(right, TreeNode):
        filepath = str(right[0])
        name = basename(filepath.split('.')[0])
        right = TreeNode(Name=name)
        right.FilePath = filepath
        right.Processed = False
        right.PollPath = None  # doesn't make sense for tips
        right.FullCommand = None
        right.EndTime = None
        right.StartTime = None
        right.TotalTime = None

    # internal node
    name = str(INTERNAL_COUNT)
    filepath = join(working_dir, name) + '.biom'
    merged = TreeNode(Name=name, Children=[left, right])
    merged.FilePath = filepath
    merged.Processed = False
    merged.PollPath = filepath + '.poll'
    merged.FullCommand = None
    merged.EndTime = None
    merged.StartTime = None
    merged.TotalTime = None

    INTERNAL_COUNT += 1
    return merged


def reset_internal_count():
    global INTERNAL_COUNT
    INTERNAL_COUNT = 0


def mergeorder(items, working_dir):
    """Code taken from http://en.literateprograms.org/Merge_sort_(Python)"""
    if len(items) < 2:
        return items
    middle = len(items) / 2
    left = mergeorder(items[:middle], working_dir)
    right = mergeorder(items[middle:], working_dir)
    return mergetree(left, right, working_dir)


def initial_nodes_to_merge(tree):
    """Determine what nodes are safe to process first

    The first nodes to process are those internal nodes that have tips as
    children
    """
    to_process = set([])
    for n in tree.tips():
        sibs_are_tips = [s.istip() for s in n.siblings()]
        if all(sibs_are_tips):
            to_process.add(n.Parent)
    return to_process


def initial_has_dependencies(tree, to_process):
    """All nodes that aren't processed up front stll have dependencies

    to_process : set of nodes that are in the first round of processing
    """
    has_dependencies = []
    for n in tree.nontips(include_self=True):
        if n not in to_process:
            has_dependencies.append(n)
    return has_dependencies


def job_complete(node, verbose=False):
    """Check if the job is complete"""
    if node.PollPath is None or node.istip():
        raise JobError("Attempting to merge tip: %s" % node.Name)

    if node.Processed:
        raise JobError("Already processed node: %s" % node.Name)

    if os.path.exists(node.PollPath):
        node.EndTime = time()
        node.TotalTime = node.EndTime - node.StartTime

        node.ExitStatus = open(node.PollPath).read().strip()
        if node.ExitStatus != '0':
            raise JobError("Node %s did not complete correctly!" % node.Name)

        if verbose:
            print "finishing %s, %f seconds" % (node.Name, node.TotalTime)

        node.Processed = True
        return True

    else:
        return False


def torque_job(cmd, pollpath, name, queue):
    """Wrap a cmd for job submission"""
    qsub_call = "qsub -k oe -N %s -q %s" % ("MOTU", queue)
    to_submit = 'echo "%s; echo $? > %s" | %s' % (cmd, pollpath, qsub_call)

    return to_submit


def local_job(cmd, pollpath, name, queue):
    """make a local job"""
    to_submit = '%s; echo $? > %s' % (cmd, pollpath)

    return to_submit


def start_job(node, merge_otus_fp, queue, wrap_call=torque_job, submit=True):
    """Starts a process"""
    strfmt = {'MergeOTUs': merge_otus_fp,
              'Output': node.FilePath,
              'BIOM_A': node.Children[0].FilePath,
              'BIOM_B': node.Children[1].FilePath}

    cmd = "%(MergeOTUs)s -i %(BIOM_A)s,%(BIOM_B)s -o %(Output)s"
    wrapped = wrap_call(cmd % strfmt, node.PollPath, node.Name, queue)

    if submit:
        system(wrapped)

    node.FullCommand = wrapped
    node.StartTime = time()
