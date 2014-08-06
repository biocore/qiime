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
        self._job_graph = nx.DiGraph()
        # self._logger: A WorkflowLogger object. Should be defined when
        # calling the subclass' _construct_job_graph method
        self._logger = None
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

    def _validate_execution_order(self, results):
        """Makes sure that the execution order represented in _job_graph has
        been respected

        Parameters
        ----------
        results : dict of {Node: AsyncResult}
            The AsyncResult objects of the executed jobs
        """
        # Adapted from
        # http://ipython.org/ipython-doc/dev/parallel/dag_dependencies.html
        self._logger.write("Validating execution order... ")
        for node in self._job_graph:
            started = results[node].metadata.started
            if started is None:
                self._logger.write("Job %s: starting time not available"
                                   % node)
                continue
            for parent in self._job_graph.predecessors(node):
                finished = results[parent].metadata.completed
                if finished is None:
                    self._logger.write("Job %s: finish time not available"
                                       % parent)
                    continue
                if started < finished:
                    self._logger.write(
                        "Job order not respected: %s should have happened "
                        "after %s\n" % (node, parent))
        self._logger.write("Done\n")

    def _validate_job_status(self, results):
        """Validates that all jobs executed finished correctly

        Parameters
        ----------
        results : dict of {Node: AsyncResult}
            The AsyncResult objects of the executed jobs
        """
        self._logger.write("\nValidating job status:\n")
        for node, ar in results.items():
            self._logger.write("\nJob %s: " % node)
            if ar.successful():
                self._logger.write("Success\n")
            else:
                self._logger.write("Error\n")
                try:
                    job_result = ar.get()
                except Exception, e:
                    job_result = e
                self._logger.write("\tJob results: %s\n"
                                   "\tPython output: %s\n"
                                   "\tStandard output: %s\n"
                                   "\tStandard error: %s\n"
                                   % (job_result, ar.pyout, ar.stdout,
                                      ar.stderr))

    def _clean_up_paths(self):
        """Removes the temporary paths"""
        if not self._retain_temp_files:
            self._logger.write("\nCleaning up temporary files")
            for fp in self._filepaths_to_remove:
                remove(fp)
            for dp in self._dirpaths_to_remove:
                rmtree(dp)

    def _job_blocker(self, results):
        # Block until all jobs are done
        self._logger.write("\nWaiting for all jobs to finish... ")
        context.wait(results.values())
        self._logger.write("Done\n")
        self._validate_job_status(results)
        self._validate_execution_order(results)
        self._clean_up_paths()
        self._logger.close()

    def __call__(self, *args, **kwargs):
        self._construct_job_graph(*args, **kwargs)

        if self._job_graph is None or self._logger is None:
            raise RuntimeError(
                "Job graph and/or logger not instantiated in the subclass")

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
            self._logger.write("Submitting job %s: %s... " % (node, job))
            if self._job_graph.node[node]['requires_deps']:
                # The job requires the results of the previous jobs, get them
                # and add them to a dict keyed by dependency node name
                deps_dict = {n: results[n].get()
                             for n in self._job_graph.predecessors(node)}
                results[node] = context.submit_async_deps(
                    deps, *job, dep_results=deps_dict)
            else:
                results[node] = context.submit_async_deps(deps, *job)
            self._logger.write("Done\n")

        if self._block:
            self._job_blocker(results)
        else:
            context.submit_async(self._job_blocker, results)
