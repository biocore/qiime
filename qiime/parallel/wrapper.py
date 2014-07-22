from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "Jai Ram Rideout",
               "Daniel McDonald", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from shutil import rmtree
from os import remove

import networkx as nx

from qiime.parallel.context import context


class ParallelWrapper(object):
    """Base class for any parallel code"""
    def __init__(self, retain_temp_files=False, block=True):
        # Set if we have to delete the temporary files or keep them
        self._retain_temp_files = retain_temp_files
        # Set if we have to wait until the job is done or we can submit
        # the waiter as another job
        self._block = block
        # self._job_graph: A networkx DAG holding the job workflow. Should be
        # defined when calling the subclass' _construct_job_graph method. Each
        # node should have two properties:
        #     - "job": a tuple with the actual job to execute in the node
        #     - "requires_deps": a boolean indicating if a dictionary with the
        #       results of the previous jobs (keyed by node name) is passed to
        #       the job using the kwargs argument with name 'dep_results'
        self._job_graph = None
        # self._log_file: The path to the log file. Should be defined when
        # calling the subclass' _construct_job_graph method
        self._log_file = None
        # Clean up variables
        self._filepaths_to_remove = []
        self._dirpaths_to_remove = []

    def _construct_job_graph(self, *args, **kwargs):
        """Constructs the workflow graph with the jobs to execute

        Raises
        ------
        NotImplementedError
            If not overwritten in a subclass
        """
        raise NotImplementedError("This method should be overwritten by the "
                                  "subclass")

    def _validate_execution_order(self, results, log_f):
        """Makes sure that the execution order represented in _job_graph has
        been respected

        Parameters
        ----------
        results : dict of {Node: AsyncResult}
            The AsyncResult objects of the executed jobs
        log_f : file like
            The open file handler.
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
