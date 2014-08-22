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

from brokit.formatdb import build_blast_db_from_fasta_path

from qiime.util import load_qiime_config
from qiime.workflow.util import generate_log_fp, WorkflowLogger
from qiime.parallel.context import context
from qiime.parallel.wrapper import ParallelWrapper
from qiime.parallel.util import (input_fasta_splitter, concatenate_files,
                                 command_wrapper, fasta_splitter_handler,
                                 blast_db_builder_handler)


class ParallelBlaster(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        """Creates the job workflow grapth to run Blast on parallel

        Parameters
        ----------
        input_fp : str
            Path to the input fasta file
        output_dir : str
            Path to the output directory. It will be created if it does
            not exists
        params : dict
            Parameters to use when calling blastall, in the form of
            {param_name: value}
        jobs_to_start : int, optional
            Number of jobs to start. Default: None - start as many jobs as
            workers in the cluster
        """
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
            job = (build_blast_db_from_fasta_path, refseqs_path, False,
                   working_dir, False)
        else:
            job = (lambda: (params['blast_db'], None), )

        self._job_graph.add_node("BUILD_BLAST_DB", job=job,
                                 requires_deps=False)
        dep_jobs.append("BUILD_BLAST_DB")

        # Split the input fasta file
        self._job_graph.add_node("SPLIT_FASTA",
                                 job=(input_fasta_splitter, input_fp,
                                      working_dir, jobs_to_start),
                                 requires_deps=False)
        dep_jobs.append("SPLIT_FASTA")

        node_names = []
        temp_outs = []
        keys = ["SPLIT_FASTA", "BUILD_BLAST_DB"]
        funcs = {"SPLIT_FASTA": fasta_splitter_handler,
                 "BUILD_BLAST_DB": blast_db_builder_handler}
        for i in range(jobs_to_start):
            node_name = "BLAST_%s" % i
            node_names.append(node_name)
            outfile = join(working_dir, "%s_out.txt" % node_name)
            temp_outs.append(outfile)
            cmd = ("%s -p blastn -m 9 -e %s -F %s -W %s -b %s %s > %s"
                   % (blastall_fp, params['e_value'], complexity_str,
                      params['word_size'], params['num_hits'], "-i %s -d %s",
                      outfile))
            self._job_graph.add_node(node_name,
                                     job=(command_wrapper, cmd, i, keys,
                                          funcs),
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
