from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import basename, splitext, abspath, join, exists
from os import makedirs
from tempfile import mkdtemp

from brokit.formatdb import build_blast_db_from_fasta_path

from qiime.align_seqs import compute_min_alignment_length
from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import (input_fasta_splitter, merge_files_from_dirs,
                                 concatenate_files, command_wrapper,
                                 fasta_splitter_handler,
                                 blast_db_builder_handler)
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp, WorkflowLogger


class ParallelAlignSeqsPyNast(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        """Creates the job workflow graph to align sequences in parallel using
        the PyNast algorithm.

        Parameters
        ----------
        input_fp : str
            Path to the input fasta file
        output_dir : str
            Path to the output directory. It will be created if it does
            not exists
        params : dict
            Parameters to use when calling align_seqs.py, in the form of
            {param_name: value}
        jobs_to_start : int, optional
            Number of jobs to start. Default: None - start as many jobs as
            workers in the cluster
        """
        # Do the parameter parsing
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        template_fp = abspath(params['template_fp'])
        blast_db = params['blast_db']
        min_length = params['min_length']

        output_dir = abspath(output_dir)

        # Create the output directory
        if not exists(output_dir):
            makedirs(output_dir)

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Generate the log file
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='align_seqs_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)
        dep_job_names = []

        if not blast_db:
            # The user did not provided the blast db, so the we have to
            # build it
            job = (build_blast_db_from_fasta_path, template_fp, False,
                   working_dir, False)
        else:
            # The user did provided the blast db, we just submit a dummy
            # job that makes the workflow easier
            job = (lambda: (blast_db, None), )

        # Add the blast db builder node to the workflow
        self._job_graph.add_node("BUILD_BLAST_DB", job=job,
                                 requires_deps=False)
        dep_job_names.append("BUILD_BLAST_DB")

        if min_length < 0:
            min_length = compute_min_alignment_length(open(input_fp, 'U'))

        # Split the input fasta file
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)
        dep_job_names.append("SPLIT_FASTA")

        # Get job commands
        output_dirs = []
        node_names = []
        keys = ["SPLIT_FASTA", "BUILD_BLAST_DB"]
        funcs = {"SPLIT_FASTA": fasta_splitter_handler,
                 "BUILD_BLAST_DB": blast_db_builder_handler}
        for i in range(jobs_to_start):
            out_dir = join(working_dir, "align_seqs_%d" % i)
            cmd = ("align_seqs.py -p %1.2f -e %d -m pynast -t %s -a %s "
                   "-o %s %s"
                   % (params['min_percent_id'], min_length, template_fp,
                      params['pairwise_alignment_method'], out_dir,
                      "-i %s -d %s"))
            output_dirs.append(out_dir)
            node_name = "AS_%d" % i
            node_names.append(node_name)
            self._job_graph.add_node(node_name,
                                     job=(command_wrapper, cmd, i, keys,
                                          funcs),
                                     requires_deps=True)
            # Adding the dependency edges to the graph
            for dep_node_name in dep_job_names:
                self._job_graph.add_edge(dep_node_name, node_name)

        # Generate the paths to the output files
        prefix = splitext(basename(input_fp))[0]
        aligned_fp = join(output_dir, "%s_aligned.fasta" % prefix)
        failures_fp = join(output_dir, "%s_failures.fasta" % prefix)
        log_fp = join(output_dir, "%s_log.txt" % prefix)

        # Merge the results by concatenating the output files
        self._job_graph.add_node("CONCAT_ALIGNED",
                                 job=(merge_files_from_dirs, aligned_fp,
                                      output_dirs, "*_aligned.fasta",
                                      concatenate_files),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_FAILURES",
                                 job=(merge_files_from_dirs, failures_fp,
                                      output_dirs, "*_failures.fasta",
                                      concatenate_files),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_LOGS",
                                 job=(merge_files_from_dirs, log_fp,
                                      output_dirs, "*_log.txt",
                                      concatenate_files),
                                 requires_deps=False)

        # Make sure that the concatenate jobs are executed after the worker
        # nodes are finished
        for node in node_names:
            self._job_graph.add_edge(node, "CONCAT_ALIGNED")
            self._job_graph.add_edge(node, "CONCAT_FAILURES")
            self._job_graph.add_edge(node, "CONCAT_LOGS")
