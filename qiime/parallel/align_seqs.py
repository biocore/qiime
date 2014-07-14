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

        # Pre-command initiation
        if not blast_db:
            # Build the blast database from reference_seqs_fp -- all procs
            # will then access one db rather than create on per proc
            blast_db, db_files_to_remove = build_blast_db_from_fasta_path(
                template_fp, output_dir=working_dir)
            self._filepaths_to_remove.extend(db_files_to_remove)

        if min_length < 0:
            min_length = compute_min_alignment_length(open(input_fp, 'U'))

        # Split the input fasta file
        fasta_fps = input_fasta_splitter(input_fp, working_dir, jobs_to_start)

        # Get job commands
        aligned_fps = []
        failures_fps = []
        log_fps = []
        node_names = []
        for i, fasta_fp in enumerate(fasta_fps):
            cmd = (
                "align_seqs.py -d %s -p %1.2f -e %d -m pynast -t %s -a %s "
                "-o %s -i %s"
                % (blast_db, params['min_percent_id'], min_length,
                   template_fp, params['pairwise_alignment_method'],
                   working_dir, fasta_fp))
            prefix = splitext(basename(fasta_fp))[0]
            aligned_fps.append(join(working_dir, "%s_aligned.fasta" % prefix))
            failures_fps.append(join(working_dir,
                                     "%s_failures.fasta" % prefix))
            log_fps.append(join(working_dir, "%s_log.txt" % prefix))
            node_name = "AS_%d" % i
            node_names.append(node_name)
            self._job_graph.add_node(node_name, job=(cmd,),
                                     requires_deps=False)

        prefix = splitext(basename(input_fp))[0]
        aligned_fp = join(output_dir, "%s_aligned.fasta" % prefix)
        failures_fp = join(output_dir, "%s_failures.fasta" % prefix)
        log_fp = join(output_dir, "%s_log.txt" % prefix)
        # Merge the results
        self._job_graph.add_node("CONCAT_ALIGNED",
                                 job=(concatenate_files, aligned_fp,
                                      aligned_fps),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_FAILURES",
                                 job=(concatenate_files, failures_fp,
                                      failures_fps),
                                 requires_deps=False)
        self._job_graph.add_node("CONCAT_LOGS",
                                 job=(concatenate_files, log_fp, log_fps),
                                 requires_deps=False)
        # Make sure that the concatenate jobs are executed after the worker
        # nodes are finished
        for node in node_names:
            self._job_graph.add_edge(node, "CONCAT_ALIGNED")
            self._job_graph.add_edge(node, "CONCAT_FAILURES")
            self._job_graph.add_edge(node, "CONCAT_LOGS")


# class ParallelAlignSeqsPyNastOld(ParallelWrapper):
#     _script_name = "align_seqs.py"
#     _job_prefix = 'ALIGN'
#     _input_splitter = ParallelWrapper._split_fasta

#     def _get_job_commands(self,
#                           fasta_fps,
#                           output_dir,
#                           params,
#                           job_prefix,
#                           working_dir,
#                           command_prefix='/bin/bash; ',
#                           command_suffix='; exit'):
#         """Generate PyNAST commands which should be submitted to cluster
#         """
#         # Create basenames for each of the output files. These will be filled
#         # in to create the full list of files created by all of the runs.
#         out_filenames = [job_prefix + '.%d_aligned.fasta',
#                          job_prefix + '.%d_failures.fasta',
#                          job_prefix + '.%d_log.txt']

#         # Initialize the command_prefix and command_suffix
#         command_prefix = command_prefix or '/bin/bash; '
#         command_suffix = command_suffix or '; exit'

#         # Create lists to store the results
#         commands = []
#         result_filepaths = []

#         # If there is a value for blast_db, pass it. If not, it
#         # will be created on-the-fly. Note that on-the-fly blast dbs
#         # are created with a string of random chars in the name, so this is safe.
#         # They shouldn't overwrite one another, and will be cleaned up.
#         if params['blast_db']:
#             blast_str = '-d %s' % params['blast_db']
#         else:
#             blast_str = ''

#         # Iterate over the input files
#         for i, fasta_fp in enumerate(fasta_fps):
#             # Each run ends with moving the output file from the tmp dir to
#             # the output_dir. Build the command to perform the move here.
#             rename_command, current_result_filepaths = self._get_rename_command(
#                 [fn % i for fn in out_filenames], working_dir, output_dir)
#             result_filepaths += current_result_filepaths

#             command = \
#                 '%s %s %s -p %1.2f -e %d -m pynast -t %s -a %s -o %s -i %s %s %s' %\
#                 (command_prefix,
#                  self._script_name,
#                  blast_str,
#                  params['min_percent_id'],
#                  params['min_length'],
#                  params['template_fp'],
#                  params['pairwise_alignment_method'],
#                  working_dir,
#                  fasta_fp,
#                  rename_command,
#                  command_suffix)

#             commands.append(command)

#         return commands, result_filepaths

#     def _write_merge_map_file(self,
#                               input_file_basename,
#                               job_result_filepaths,
#                               params,
#                               output_dir,
#                               merge_map_filepath):

#         f = open(merge_map_filepath, 'w')

#         out_filepaths = [
#             '%s/%s_aligned.fasta' % (output_dir, input_file_basename),
#             '%s/%s_failures.fasta' % (output_dir,
#                                       input_file_basename),
#             '%s/%s_log.txt' % (output_dir, input_file_basename)]

#         aligned_fps = []
#         failures_fps = []
#         log_fps = []

#         for fp in job_result_filepaths:
#             if fp.endswith('_aligned.fasta'):
#                 aligned_fps.append(fp)
#             elif fp.endswith('_failures.fasta'):
#                 failures_fps.append(fp)
#             else:
#                 log_fps.append(fp)

#         for in_files, out_file in\
#                 zip([aligned_fps, failures_fps, log_fps], out_filepaths):
#             f.write('\t'.join(in_files + [out_file]))
#             f.write('\n')
#         f.close()
