from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import basename, splitext, abspath, join
from tempfile import mkdtemp

import networkx as nx
from brokit.formatdb import build_blast_db_from_fasta_path

from qiime.align_seqs import compute_min_alignment_length
from qiime.parallel.util import (ParallelWrapper, input_fasta_splitter,
                                 concatenate_files)
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp
from qiime.util import get_qiime_temp_dir


def command_wrapper(cmd, idx, dep_results=None):
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

    if "BUILD_BLAST_DB" not in dep_results:
        raise ValueError("Wrong job graph workflow. Node 'BUILD_BLAST_DB' "
                         "not listed as dependency of current node")
    blast_db, db_files_to_remove = dep_results["BUILD_BLAST_DB"]
    cmd = cmd % (fasta_fps[idx], blast_db)
    return system_call(cmd)


def concatenate_wrapper(output_fp, output_dirs, suffix):
    """Wraps the concatenate_files function so it can generate the actual list
    of files to concatenate

    Parameters
    ----------
    output_fp : str
        The path to the output file
    output_dirs : list of str
        The list of directories in which we should search for the files
    suffix : str
        The suffix of the files that we have to search for. It should already
        include the '*' if needed
    """
    # Importing glob here so it is available to the workers
    from glob import glob
    files = []
    for out_dir in output_dirs:
        files.extend(glob(join(out_dir, suffix)))
    concatenate_files(output_fp, files)


class ParallelAlignSeqsPyNast(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()
        # Do the parameter parsing
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        template_fp = abspath(params['template_fp'])
        blast_db = params['blast_db']
        min_length = params['min_length']

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Generate the log file
        self._log_file = generate_log_fp(output_dir)

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='align_seqs_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)
        dep_job_names = []

        # Pre-command initiation
        if not blast_db:
            # Build the blast database from reference_seqs_fp -- all procs
            # will then access one db rather than create on per proc
            self._job_graph.add_node("BUILD_BLAST_DB",
                                     job=(build_blast_db_from_fasta_path,
                                          template_fp, False, working_dir,
                                          False),
                                     requires_deps=False)
        else:
            func = lambda: (blast_db, None)
            self._job_graph.add_node("BUILD_BLAST_DB",
                                     job=(func, ),
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
        for i in range(jobs_to_start):
            # out_dir = mkdtemp(prefix="align_seqs_%d" % i, dir=wor)
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
                                     job=(command_wrapper, cmd, i),
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
                                 job=(concatenate_wrapper, aligned_fp,
                                      output_dirs, "*_aligned.fasta"),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_FAILURES",
                                 job=(concatenate_wrapper, failures_fp,
                                      output_dirs, "*_failures.fasta"),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_LOGS",
                                 job=(concatenate_wrapper, log_fp, output_dirs,
                                      "*_log.txt"),
                                 requires_deps=False)
        # Make sure that the concatenate jobs are executed after the worker
        # nodes are finished
        for node in node_names:
            self._job_graph.add_edge(node, "CONCAT_ALIGNED")
            self._job_graph.add_edge(node, "CONCAT_FAILURES")
            self._job_graph.add_edge(node, "CONCAT_LOGS")
