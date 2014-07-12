#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import split, abspath, basename, splitext, join
from shutil import copyfile
from tempfile import mkdtemp

import networkx as nx
from brokit.formatdb import build_blast_db_from_fasta_path
from skbio.parse.sequences import parse_fasta

from qiime.identify_chimeric_seqs import make_cidx_file
from qiime.parse import parse_tmp_to_final_filepath_map_file
from qiime.util import write_degapped_fasta_to_file, count_seqs
from qiime.parallel.util import ParallelWrapper, qiime_config
from qiime.parallel.context import context
from qiime.split import split_fasta
from qiime.workflow.util import generate_log_fp


def _concatenate_output(output_fp, temp_out_fps):
    with open(output_fp, 'w') as out_f:
            for tmp_fp in temp_out_fps:
                with open(tmp_fp, 'U') as in_f:
                    for line in in_f:
                        out_f.write(line)


def _split_fasta(input_fp, output_dir, num):
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


class ParallelChimericSequenceIdentifier(ParallelWrapper):

    def _construct_job_graph(self, input_fp, output_dir, params,
                             jobs_to_start=None):
        self._job_graph = nx.DiGraph()

        # Do the parameter parsing
        input_fp = abspath(input_fp)
        output_dir = abspath(output_dir)
        if jobs_to_start is None:
            # default to the number of workers
            jobs_to_start = context.get_number_of_workers()
        method = params['chimera_detection_method']
        self._log_file = generate_log_fp(output_dir)

        # Get a temp folder
        working_dir = mkdtemp(prefix='CSI_', dir=output_dir)
        self._dirpaths_to_remove.append(working_dir)

        if method == 'blast_fragments':
            ref_seqs_fp = abspath(params['reference_seqs_fp'])
            # We first need to create the blast database
            handler = context.submit_async(build_blast_db_from_fasta_path,
                                           ref_seqs_fp, output_dir=working_dir)
            blast_db, db_files_to_remove = handler.get()
            self._filepaths_to_remove.extend(db_files_to_remove)
            params['blast_db'] = blast_db
            # Split the input_fp in multiple files
            handler = context.submit_async(_split_fasta, input_fp, working_dir,
                                           jobs_to_start)
            fasta_fps = handler.get()
            self._filepaths_to_remove.extend(fasta_fps)
            # Now we create all the commands
            temp_out_fps = []
            ics_nodes = []
            for i, fasta_fp in enumerate(fasta_fps):
                # Create the output path
                temp_out_fp = join(
                    working_dir,
                    "%s_%d.txt" % (splitext(basename(fasta_fp))[0], i))
                temp_out_fps.append(temp_out_fp)
                cmd = (
                    "identify_chimeric_seqs.py -i %s -t %s -m blast_fragments "
                    "-o %s -n %s -d %s -e %s -b %s"
                    % (fasta_fp, params['id_to_taxonomy_fp'], temp_out_fp,
                       params['num_fragments'], params['taxonomy_depth'],
                       params['max_e_value'], params['blast_db']))
                node_name = "ICS_%d" % i
                ics_nodes.append(node_name)
                self._job_graph.add_node(node_name, job=(cmd,))

            # Adding temp_out_fps to the clean-up variable
            self._filepaths_to_remove.extend(temp_out_fps)
            # Merge the results by concatenating the output files
            output_fp = join(output_dir, "chimeric_seqs.txt")
            self._job_graph.add_node("CONCAT", job=(_concatenate_output,
                                                    output_fp, temp_out_fps))
            # Make sure that the "CONCAT" node is the latest executed node
            for node in ics_nodes:
                self._job_graph.add_edge(node, "CONCAT")

        elif method == 'ChimeraSlayer':
            # We first need to copy the reference files to the working
            # directory because ChimeraSlayer creates an index file of the
            # reference and it will crash without write permission in the
            # reference sequences directory
            aln_ref_seqs_fp = abspath(params['aligned_reference_seqs_fp'])
            dest_aln_ref_seqs_fp = join(working_dir, basename(aln_ref_seqs_fp))
            copyfile(aln_ref_seqs_fp, dest_aln_ref_seqs_fp)
            aln_ref_seqs_fp = dest_aln_ref_seqs_fp
            # Adding the copy to the clean up variable
            self._filepaths_to_remove.append(aln_ref_seqs_fp)

            ref_seqs_fp = params['reference_seqs_fp']
            if ref_seqs_fp:
                ref_seqs_fp = abspath(ref_seqs_fp)
                # Copy the unaligned reference database
                dest_ref_seqs_fp = join(working_dir, basename(ref_seqs_fp))
                copyfile(ref_seqs_fp, dest_ref_seqs_fp)
                ref_seqs_fp = dest_ref_seqs_fp
            else:
                # The unaligned reference database does not exists, create it
                with open(aln_ref_seqs_fp, 'U') as f:
                    ref_seqs_fp = write_degapped_fasta_to_file(
                        parse_fasta(f), tmp_dir=working_dir)

            # Add the copied/new ref_seqs_fp to the clean up variable
            self._filepaths_to_remove.append(ref_seqs_fp)

            # Build the blast db of reference, otherwise ChimeraSlayer will do
            # it and parallel jobs clash
            handler = context.submit_async(build_blast_db_from_fasta_path,
                                           ref_seqs_fp, output_dir=working_dir)
            blast_db, db_files_to_remove = handler.get()
            self._filepaths_to_remove.extend(db_files_to_remove)

            # Make the index file. ChimeraSlayer first checks to see if the
            # index file exists. If not, it tries to create it. This can lead
            # to race conditions if several parallel jobs try to create it at
            # the same time
            self._job_graph.add_node("MAKE_CIDX", job=(make_cidx_file,
                                                       aln_ref_seqs_fp))
            self._filepaths_to_remove.append("%s.cidx" % aln_ref_seqs_fp)

            # Split the input_fp in multiple files
            handler = context.submit_async(_split_fasta, input_fp, working_dir,
                                           jobs_to_start)
            fasta_fps = handler.get()
            self._filepaths_to_remove.extend(fasta_fps)

            min_div_ratio_str = ("--min_div_ratio %s" % params['min_div_ratio']
                                 if params['min_div_ratio'] else "")

            # Create the commands
            temp_out_fps = []
            ics_nodes = []
            for i, fasta_fp in enumerate(fasta_fps):
                # Create the output path
                temp_out_fp = join(
                    working_dir,
                    "%s_%d.txt" % (splitext(basename(fasta_fp))[0], i))
                temp_out_fps.append(temp_out_fp)
                cmd = (
                    "identify_chimeric_seqs.py -i %s -a %s -m ChimeraSlayer "
                    "-o %s -r %s %s"
                    % (fasta_fp, aln_ref_seqs_fp, temp_out_fp,
                       ref_seqs_fp, min_div_ratio_str))
                node_name = "ICS_%d" % i
                ics_nodes.append(node_name)
                self._job_graph.add_node(node_name, job=(cmd,))
                # This should be executed after the "MAKE_CIDX" node
                self._job_graph.add_edge("MAKE_CIDX", node_name)

            # Adding temp_out_fps to the clean_up variable
            self._filepaths_to_remove.extend(temp_out_fps)
            # Merge the results by concatenating the output files
            output_fp = join(output_dir, "chimeric_seqs.txt")
            self._job_graph.add_node("CONCAT", job=(_concatenate_output,
                                                    output_fp, temp_out_fps))
            # Make sure that the "CONCAT" node is the latest executed node
            for node in ics_nodes:
                self._job_graph.add_edge(node, "CONCAT")

        else:
            raise ValueError(
                "Unrecognized chimera detection method '%s'." % method)


# class ParallelChimericSequenceIdentifierOld(ParallelWrapper):
#     _script_name = 'identify_chimeric_seqs.py'
#     _input_splitter = ParallelWrapper._split_fasta
#     _job_prefix = 'CHIM'
#     _process_run_results_f = \
#         'qiime.parallel.identify_chimeric_seqs.basic_process_run_results_f'

#     def _precommand_initiation(self, input_fp, output_dir, working_dir,
#                                params):
#         if params['chimera_detection_method'] == 'blast_fragments':
#             blast_db, db_files_to_remove = \
#                 build_blast_db_from_fasta_path(params['reference_seqs_fp'],
#                                                output_dir=working_dir)
#             self.files_to_remove += db_files_to_remove
#             params['blast_db'] = blast_db
#         elif params['chimera_detection_method'] == 'ChimeraSlayer':
#             # copy the reference files to working dir
#             # ChimeraSlayer creates an index file of the ref and
#             # will crash without write permission in the ref seqs dir
#             aligned_reference_seqs_fp = params['aligned_reference_seqs_fp']
#             _, new_ref_filename = split(aligned_reference_seqs_fp)
#             copy(aligned_reference_seqs_fp, working_dir)
#             aligned_reference_seqs_fp = working_dir + "/" + new_ref_filename

#             self.files_to_remove.append(aligned_reference_seqs_fp)
#             params['aligned_reference_seqs_fp'] = aligned_reference_seqs_fp

#             # if given, also copy the unaligned ref db
#             reference_seqs_fp = params['reference_seqs_fp']
#             if reference_seqs_fp:
#                 _, new_ref_filename = split(reference_seqs_fp)
#                 copy(reference_seqs_fp, working_dir)
#                 reference_seqs_fp = working_dir + "/" + new_ref_filename
#             else:
#                 # otherwise create it
#                 reference_seqs_fp = write_degapped_fasta_to_file(
#                     parse_fasta(open(aligned_reference_seqs_fp)),
#                     tmp_dir=working_dir)
#             # delete it afterwards
#             self.files_to_remove.append(reference_seqs_fp)
#             params['reference_seqs_fp'] = reference_seqs_fp

#             # build blast db of reference, otherwise ChimeraSlayer will do it
#             # and parallel jobs clash
#             _, db_files_to_remove = \
#                 build_blast_db_from_fasta_path(reference_seqs_fp)
#             self.files_to_remove += db_files_to_remove

#             # make the index file globally
#             # Reason: ChimeraSlayer first checks to see if the index file is
#             # there. If not it tries to create it. This can lead to race
#             # condition if several parallel jobs try to create it at the same
#             # time.
#             make_cidx_file(aligned_reference_seqs_fp)
#             self.files_to_remove.append(aligned_reference_seqs_fp + ".cidx")
#         else:
#             raise ValueError("Unrecognized chimera detection method '%s'." %
#                              params['chimera_detection_method'])

#     def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
#                           working_dir, command_prefix='/bin/bash; ',
#                           command_suffix='; exit'):
#         """Generate identify_chimeric_seqs.py commands which should be run."""
#         # Create basenames for each of the output files. These will be filled
#         # in to create the full list of files created by all of the runs.
#         out_filenames = [job_prefix + '.%d_chimeric.txt']

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

#             optional_options = ""
#             if params['chimera_detection_method'] == 'blast_fragments':
#                 command = \
#                     '%s %s -i %s -t %s -m blast_fragments -o %s -n %s -d %s -e %s -b %s %s %s' % \
#                     (command_prefix,
#                      self._script_name,
#                      fasta_fp,
#                      params['id_to_taxonomy_fp'],
#                      working_dir + "/" + out_filenames[0] % i,
#                      params['num_fragments'],
#                      params['taxonomy_depth'],
#                      params['max_e_value'],
#                      params['blast_db'],
#                      rename_command,
#                      command_suffix)
#             elif params['chimera_detection_method'] == 'ChimeraSlayer':
#                 optional_options = ""
#                 if params['min_div_ratio']:
#                     optional_options += " --min_div_ratio %s" % \
#                                         params['min_div_ratio']
#                 if params['reference_seqs_fp']:
#                     optional_options += " -r %s" % params['reference_seqs_fp']

#                 command = \
#                     '%s %s -i %s -a %s -m ChimeraSlayer -o %s %s %s %s' % \
#                     (command_prefix,
#                      self._script_name,
#                      fasta_fp,
#                      params['aligned_reference_seqs_fp'],
#                      working_dir + "/" + out_filenames[0] % i,
#                      optional_options,
#                      rename_command,
#                      command_suffix)
#             else:
#                 raise NotImplementedError
#             commands.append(command)
#         return commands, result_filepaths

#     def _get_poller_command(self,
#                             expected_files_filepath,
#                             merge_map_filepath,
#                             deletion_list_filepath,
#                             command_prefix='/bin/bash; ',
#                             command_suffix='; exit'):
#         """Generate command to initiate a poller to monitior/process completed runs
#         """
#         result = '%s poller.py -f %s -p %s -m %s -d %s -t %d %s' % \
#             (command_prefix,
#              expected_files_filepath,
#              self._process_run_results_f,
#              merge_map_filepath,
#              deletion_list_filepath,
#              self._seconds_to_sleep,
#              command_suffix)
#         return result, []

#     def _write_merge_map_file(self, input_file_basename, job_result_filepaths,
#                               params, output_dir, merge_map_filepath,
#                               failures=False):
#         f = open(merge_map_filepath, 'w')
#         out_filepaths = [params['output_fp']]

#         chims_fps = []
#         logs_fps = []  # logs_fp currently not used

#         for fp in job_result_filepaths:
#             if fp.endswith('_chimeric.txt'):
#                 chims_fps.append(fp)
#             else:
#                 logs_fps.append(fp)

#         for in_files, out_file in zip([chims_fps], out_filepaths):
#             f.write('\t'.join(in_files + [out_file]))
#             f.write('\n')
#         f.close()


# def basic_process_run_results_f(f):
#     """ Copy each list of infiles to each outfile and delete infiles

#         f: file containing one set of mapping instructions per line

#         example f:
#          f1.txt f2.txt f3.txt f_combined.txt
#          f1.log f2.log f3.log f_combined.log

#         If f contained the two lines above, this function would
#          concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
#          and f1.log, f2.log, and f3.log into f_combined.log
#     """
#     infiles_lists, out_filepaths = parse_tmp_to_final_filepath_map_file(f)
#     for infiles_list, out_filepath in zip(infiles_lists, out_filepaths):
#         try:
#             of = open(out_filepath, 'w')
#         except IOError:
#             raise IOError("Poller can't open final output file: %s" % out_filepath +
#                           "\nLeaving individual jobs output.\n Do you have write access?")

#         for fp in infiles_list:
#             for line in open(fp):
#                 of.write('%s\n' % line.strip('\n'))
#         of.close()
#     # It is a good idea to have your clean_up_callback return True.
#     # That way, if you get mixed up and pass it as check_run_complete_callback,
#     # you'll get an error right away rather than going into an infinite loop
#     return True
