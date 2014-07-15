#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import abspath, exists, basename, splitext
from os import makedirs
from tempfile import mkdtemp

from brokit.formatdb import build_blast_db_from_fasta_path

from qiime.parallel.util import ParallelWrapper, input_fasta_splitter
from qiime.parallel.context import context
from qiime.workflow.util import generate_log_fp


class ParallelTaxonomyAssigner(ParallelWrapper):
    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        # Create the workflow graph
        self._job_graph = nx.DiGraph()

        # Create the output directory if it does not exists
        output_dir = abspath(output_dir)
        if not exists(output_dir):
            makedirs(output_dir)

        # If the number of jobs to start is not provided, we default to the
        # number of workers
        if jobs_to_start is None:
            jobs_to_start = context.get_number_of_workers()

        # Generate the log file
        self._log_file = genreate_log_fp(output_dir)

        # Get a folder to store the temporary files
        working_dir = mkdtemp(prefix='tax_assigner_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        # PRECOMMAND - ONLY BLAST PARALLEL ASSIGNER NEEDS IT
        # TODO

        # Split the input fasta file
        # self._job_graph.add_node("SPLIT_FASTA",
        #                          job=(input_fasta_splitter, input_fp,
        #                               working_dir, jobs_to_start),
        #                          requires_deps=False)
        fasta_fps = input_fasta_splitter(input_fp, working_dir, jobs_to_start)

        # Build the commands
        output_dirs = []
        node_names = []
        for cmd, node_name, out_dir in self._cmd_generator(fasta_fps,
                                                           working_dir,
                                                           params):
            output_dirs.append(out_dir)
            node_names.append(node_name)
            self._job_graph.add_node(node_name, job=(cmd, ),
                                     requires_deps=False)

        # Merge the results
        prefix = splitext(basename(input_fp))[0]
        out_tax_fp = join(output_dir, "%s_tax_assignments.txt" % prefix)
        log_fp = join(output_dir, "%s_tax_assignments.log" % prefix)


class ParallelBlastTaxonomyAssigner(ParallelTaxonomyAssigner):
    def _cmd_generator(self, fasta_fps, working_dir, params):
        for i, fasta_fp in enumerate(fasta_fps):
            node_name = "PTA_%d" % i
            output_dir = join(working_dir, node_name)
            cmd = (
                "assign_taxonomy.py -o %s -i %s -m blast -e %s -b %s -t %s"
                % (output_dir, fasta_fp, params['e_value'],
                   params['blast_db'], params['id_to_taxonomy_fp']))
            yield cmd, node_name, output_dir


# class ParallelTaxonomyAssigner(ParallelWrapper):
#     _script_name = 'assign_taxonomy.py'
#     _input_splitter = ParallelWrapper._split_fasta

#     def _build_job_commands(self, tax_specific_param_str, fasta_fps,
#                             output_dir, params, job_prefix, working_dir,
#                             command_prefix='/bin/bash; ',
#                             command_suffix='; exit'):
#         """Generate assign_taxonomy.py commands which should be run."""
#         # Create basenames for each of the output files. These will be filled
#         # in to create the full list of files created by all of the runs.
#         out_filenames = [job_prefix + '.%d_tax_assignments.log',
#                          job_prefix + '.%d_tax_assignments.txt']

#         # Create lists to store the results.
#         commands = []
#         result_filepaths = []

#         # Iterate over the input files.
#         for i, fasta_fp in enumerate(fasta_fps):
#             # Each run ends with moving the output file from the tmp dir to
#             # the output_dir. Build the command to perform the move here.
#             rename_command, current_result_filepaths = \
#                 self._get_rename_command([fn % i for fn in out_filenames],
#                                          working_dir, output_dir)
#             result_filepaths += current_result_filepaths

#             command = '%s %s %s -o %s -i %s %s %s' %\
#                 (command_prefix,
#                  self._script_name,
#                  tax_specific_param_str,
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
#                               merge_map_filepath,
#                               failures=False):
#         """
#         """
#         f = open(merge_map_filepath, 'w')

#         out_filepaths = [
#             '%s/%s_tax_assignments.txt' % (output_dir, input_file_basename),
#             '%s/%s_tax_assignments.log' % (output_dir, input_file_basename)]

#         assignment_fps = []
#         log_fps = []

#         for fp in job_result_filepaths:
#             if fp.endswith('_tax_assignments.txt'):
#                 assignment_fps.append(fp)
#             else:
#                 log_fps.append(fp)

#         for in_files, out_file in\
#                 zip([assignment_fps, log_fps], out_filepaths):
#             f.write('\t'.join(in_files + [out_file]))
#             f.write('\n')
#         f.close()


# class ParallelRdpTaxonomyAssigner(ParallelTaxonomyAssigner):
#     _job_prefix = 'RDP'

#     def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
#                           working_dir, command_prefix=None,
#                           command_suffix='; exit'):
#         command_prefix = command_prefix or \
#             '/bin/bash; export RDP_JAR_PATH=%s; ' % params['rdp_classifier_fp']

#         rdp_params = '-m rdp -c %1.2f --rdp_max_memory %d ' % (
#             params['confidence'], params['rdp_max_memory'])
#         if params['id_to_taxonomy_fp'] and params['reference_seqs_fp']:
#             rdp_params += '-t %s -r %s' % (params['id_to_taxonomy_fp'],
#                                            params['reference_seqs_fp'])

#         return self._build_job_commands(rdp_params, fasta_fps, output_dir,
#                                         params, job_prefix, working_dir, command_prefix,
#                                         command_suffix)


# class ParallelUclustConsensusTaxonomyAssigner(ParallelTaxonomyAssigner):
#     _job_prefix = 'UCTA'

#     def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
#                           working_dir, command_prefix=None,
#                           command_suffix='; exit'):
#         command_prefix = command_prefix or ''

#         uclust_params = ' '.join(
#             ['-m uclust',
#              '--uclust_min_consensus_fraction %f' % params['uclust_min_consensus_fraction'],
#              '--uclust_similarity %f' % params['uclust_similarity'],
#              '--uclust_max_accepts %d' % params['uclust_max_accepts'],
#              '-t %s' % params['id_to_taxonomy_fp'],
#              '-r %s' % params['reference_seqs_fp']])

#         return self._build_job_commands(uclust_params,
#                                         fasta_fps,
#                                         output_dir,
#                                         params,
#                                         job_prefix,
#                                         working_dir,
#                                         command_prefix,
#                                         command_suffix)


# class ParallelBlastTaxonomyAssigner(ParallelTaxonomyAssigner):
#     _job_prefix = 'BTA'

#     def _precommand_initiation(
#             self, input_fp, output_dir, working_dir, params):
#         if not params['blast_db']:
#             # Build the blast database from the reference_seqs_fp -- all procs
#             # will then access one db rather than create one per proc.
#             blast_db, db_files_to_remove = \
#                 build_blast_db_from_fasta_path(params['reference_seqs_fp'])
#             self.files_to_remove += db_files_to_remove
#             params['blast_db'] = blast_db

#     def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
#                           working_dir, command_prefix=None,
#                           command_suffix='; exit'):
#         command_prefix = command_prefix or \
#             '/bin/bash; cd %s; export BLASTMAT=%s;' % (working_dir,
#                                                        params['blastmat_dir'])

#         blast_params = '-m blast -e %s -b %s -t %s ' % (
#             params['e_value'], params['blast_db'],
#             params['id_to_taxonomy_fp'])

#         return self._build_job_commands(blast_params, fasta_fps, output_dir,
#                                         params, job_prefix, working_dir, command_prefix,
#                                         command_suffix)
