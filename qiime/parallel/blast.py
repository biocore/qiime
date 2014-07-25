#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import splitext, abspath, basename, exists, join
from os import makedirs
from tempfile import mkdtemp

import networkx as nx
from brokit.formatdb import build_blast_db_from_fasta_path

from qiime.util import load_qiime_config
from qiime.workflow.util import generate_log_fp, WorkflowLogger
from qiime.parallel.context import context
from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import input_fasta_splitter, concatenate_files


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


class ParallelBlaster(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        # Do the parameter parsing
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        complexity_str = ('T' if not params['disable_low_complexity_filter']
                          else 'F')
        output_dir = abspath(output_dir)
        # Create the output directory
        if not exists(output_dir):
            makedirs(output_dir)

        # Generate the log file
        self._logger = WorkflowLogger(generate_log_fp(output_dir))

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Get the blastall executable form the qiime config
        blastall_fp = load_qiime_config()['blastall_fp']

        # Get a working directory to store the temporary files
        working_dir = mkdtemp(prefix='blast_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        dep_jobs = []
        if params['refseqs_path']:
            # Build the blast database from the refseqs_path -- all procs
            # will then access one db rather than create one per proc
            refseqs_path = abspath(params['refseqs_path'])
            self._job_graph.add_node("BUILD_BLAST_DB",
                                     job=(build_blast_db_from_fasta_path,
                                          refseqs_path, False, working_dir,
                                          False),
                                     requires_deps=False)
        else:
            func = lambda: (params['blast_db'], None)
            self._job_graph.add_node("BUILD_BLAST_DB", job=(func, ),
                                     requires_deps=False)
            # blast_db, db_files_to_remove = build_blast_db_from_fasta_path(
            #     refseqs_path, False, working_dir, False)
            # params['blast_db'] = blast_db
        dep_jobs.append("BUILD_BLAST_DB")

        # Split the input fasta file
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)
        dep_jobs.append("SPLIT_FASTA")

        node_names = []
        temp_outs = []
        for i in range(jobs_to_start):
            node_name = "BLAST_%s" % i
            node_names.append(node_name)
            outfile = join(working_dir, "%s_out.txt" % node_name)
            temp_outs.append(outfile)
            cmd = ("%s -p blastn -m 9 -e %s -F %s -W %s -b %s %s > %s"
                   % (blastall_fp, params['e_value'], complexity_str,
                      params['word_size'], params['num_hits'], "-i %s -d %s",
                      outfile))
            self._job_graph.add_node(node_name, job=(command_wrapper, cmd, i),
                                     requires_deps=True)
            for job in dep_jobs:
                self._job_graph.add_edge(job, node_name)

        # Merge results by concatenating output files
        prefix = splitext(basename(input_fp))[0]
        output_fp = join(output_dir, "%s_blast_out.txt" % prefix)
        self._job_graph.add_node("CONCAT",
                                 job=(concatenate_files, output_fp, temp_outs),
                                 requires_deps=False)
        for node_name in node_names:
            self._job_graph.add_edge(node_name, "CONCAT")
