#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "Jai Ram Rideout",
               "Daniel McDonald", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import remove
from shutil import rmtree

import networkx as nx

from qiime.util import load_qiime_config
from qiime.parallel.context import context


qiime_config = load_qiime_config()


class ParallelWrapper(object):
    """"""
    def __init__(self, retain_temp_files=False, block=True):
        self._retain_temp_files = retain_temp_files
        self._block = block
        # These attributes should be defined when calling the subclass'
        # _construct_job_graph method
        # A networkx DAG holding the job workflow. Each node should have two
        # properties "job" and "requires_deps". Job holds a tuple with the
        # actual job to execute in the node. Requires_deps is a boolean. If
        # true, a dictionary with the results of the previous jobs (keyed by
        # node name) is passed to the job using the kwargs argument with
        # name dep_results
        self._job_graph = None
        self._log_file = None
        # Clean up variables
        self._filepaths_to_remove = []
        self._dirpaths_to_remove = []

    def _construct_job_graph(self, **kwargs):
        """Constructs the workflow graph and the jobs to execute"""
        raise NotImplementedError(
            "This method should be overwritten by the subclass")

    def _validate_execution_order(self, results, log_f):
        """Makes sure that the execution order represented in _job_graph has
        been respected

        Parameters
        ----------
        results : dict of {Node: AsyncResult}
            The AsyncResult objects of the executed jobs
        """
        # Adapted from
        # http://ipython.org/ipython-doc/dev/parallel/dag_dependencies.html
        log_f.write("Validating execution order... ")
        for node in self._job_graph:
            started = results[node].metadata.started
            if started is None:
                log_f.write("Job %s: metadata not available" % node)
                continue
            for parent in self._job_graph.predecessors(node):
                finished = results[parent].metadata.completed
                if finished is None:
                    log_f.write("Job %s: metadata not available" % parent)
                    continue
                if started < finished:
                    log_f.write(
                        "Job order not respected: %s should have happened "
                        "after %s\n" % (node, parent))
        log_f.write("Done\n")

    def _validate_job_status(self, results, log_f):
        """Validates that all jobs executed finished correctly

        Parameters
        ----------
        results : dict of {Node: AsyncResult}
            The AsyncResult objects of the executed jobs
        log_f : file object
            The open log file handler
        """
        log_f.write("\nValidating job status:\n")
        for node, ar in results.items():
            log_f.write("\nJob %s: " % node)
            if ar.successful():
                log_f.write("Success\n")
            else:
                log_f.write("Error\n")
                try:
                    job_result = ar.get()
                except Exception, e:
                    job_result = e
                log_f.write("\tJob results: %s\n"
                            "\tPython output: %s\n"
                            "\tStandard output: %s\n"
                            "\tStandard error: %s\n"
                            % (job_result, ar.pyout, ar.stdout, ar.stderr))

    def _clean_up_paths(self):
        """Removes the temporary paths"""
        if not self._retain_temp_files:
            for fp in self._filepaths_to_remove:
                remove(fp)
            for dp in self._dirpaths_to_remove:
                rmtree(fp)

    def _job_blocker(self, results, log_f):
        # Block until all jobs are done
        log_f.write("\nWaiting for all jobs to finish... ")
        context.wait(results.values())
        log_f.write("Done\n")
        self._validate_job_status(results, log_f)
        self._validate_execution_order(results, log_f)
        log_f.close()

    def __call__(self, *args, **kwargs):
        self._construct_job_graph(*args, **kwargs)

        if self._job_graph is None or self._log_file is None:
            raise RuntimeError(
                "Job graph and log file not instantiated in the subclass")

        log_f = open(self._log_file, 'w')

        # We need to submit the jobs to ipython in topological order, so we can
        # actually define the dependencies between jobs. Adapted from
        # http://ipython.org/ipython-doc/dev/parallel/dag_dependencies.html
        results = {}
        for node in nx.topological_sort(self._job_graph):
            # Get the list of predecessor jobs
            deps = [results[n] for n in self._job_graph.predecessors(node)]
            # Get the tuple with the job to run
            job = self._job_graph.node[node]['job']
            # We can now submit the job taking into account the dependencies
            log_f.write("Submitting job %s: %s... " % (node, job))
            if self._job_graph.node[node]['requires_deps']:
                # The job requires the results of the previous jobs, get them
                # and add them to a dict keyed by dependency node name
                deps_dict = {n: results[n].get()
                             for n in self._job_graph.predecessors(node)}
                results[node] = context.submit_async_deps(
                    deps, *job, dep_results=deps_dict)
            else:
                results[node] = context.submit_async_deps(deps, *job)
            log_f.write("Done\n")

        if self._block:
            self._job_blocker(results, log_f)
        else:
            context.submit_async(self._job_blocker, results, log_f)


def concatenate_files(output_fp, temp_out_fps):
    with open(output_fp, 'w') as out_f:
        for tmp_fp in temp_out_fps:
            with open(tmp_fp, 'U') as in_f:
                for line in in_f:
                    out_f.write(line)


def merge_files_from_dirs(output_fp, output_dirs, format_str, merge_func):
    """Wraps the concatenate_files function so it can generate the actual list
    of files to concatenate from the directories in output_dirs

    Parameters
    ----------
    output_fp : str
        The path to the output file
    output_dirs : list of str
        The list of directories in which we should search for the files
    format_str : str
        The formatted string of the files that we have to search for. It should
        include any wildcard that glob can parse (e.g. '*')
    merge_func : function
        The function used to merge the results. Signature: f(output_fp, files)
    """
    # Importing here so it is available to the workers
    from glob import glob
    from os.path import join
    files = []
    for out_dir in output_dirs:
        files.extend(glob(join(out_dir, format_str)))
    merge_func(output_fp, files)


def input_fasta_splitter(input_fp, output_dir, num):
    # Importing here so it becomes available on the workers
    from os.path import basename, splitext
    from qiime.util import count_seqs
    from qiime.split import split_fasta
    # First compute the number of sequences per file
    # Count the number of sequences in the fasta file
    num_input_seqs = count_seqs(input_fp)[0]

    # divide the number of sequences by the number of jobs to start
    num_seqs_per_file = num_input_seqs / num

    # if we don't have a perfect split, round up
    if num_seqs_per_file % 1 != 0:
        num_seqs_per_file += 1

    # Get the number of sequences as an integer
    num_seqs_per_file = int(num_seqs_per_file)

    # Generate a prefix for the files
    prefix = splitext(basename(input_fp))[0]
    fasta_fps = split_fasta(open(input_fp), num_seqs_per_file, prefix,
                            working_dir=output_dir)
    return fasta_fps


class BufferedWriter():

    """A file like object that delays writing to file without keeping an open filehandle

    This class comes useful in scenarios were potentially many open fhs are needed
    (e.g. during splitting of inputs for parallelization). Since
    each OS limits the max number of open fh at any time, we provide a fh like class that
    can be used much like a regular (writable) fh, but without keeping the fh open permanently.
    Using a larger buffer size speeds up the program by using less of the expensive open/close
    IO operations.
    """

    def __init__(self, filename, buf_size=100):
        """
        filename: name of file to write to in append mode

        buf_size: buffer size in chunks. Each write operations counts as one chunk.
        """

        if(buf_size < 1):
            raise ValueError("Invalid buf_size. Must be 1 or larger.")

        self.buffer = []
        self.buf_size = buf_size
        self.filename = filename

        # touch the file
        fh = open(self.filename, "w")
        fh.close()

    def __del__(self):
        self._flush()

    def close(self):
        self._flush()

    def write(self, line):
        """write line to BufferedWriter"""

        self.buffer.append(line)
        if (len(self.buffer) > self.buf_size):
            self._flush()

    def _flush(self):
        """Write buffer to file"""

        fh = open(self.filename, "a")
        fh.write("".join(self.buffer))
        fh.close()

        self.buffer = []
