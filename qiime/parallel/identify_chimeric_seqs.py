#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import split
from shutil import copy

from bfillings.formatdb import build_blast_db_from_fasta_path

from skbio.parse.sequences import parse_fasta

from qiime.identify_chimeric_seqs import make_cidx_file
from qiime.parse import parse_tmp_to_final_filepath_map_file
from qiime.util import write_degapped_fasta_to_file
from qiime.parallel.util import ParallelWrapper


class ParallelChimericSequenceIdentifier(ParallelWrapper):
    _script_name = 'identify_chimeric_seqs.py'
    _input_splitter = ParallelWrapper._split_fasta
    _job_prefix = 'CHIM'
    _process_run_results_f = \
        'qiime.parallel.identify_chimeric_seqs.basic_process_run_results_f'

    def _precommand_initiation(self, input_fp, output_dir, working_dir,
                               params):
        if params['chimera_detection_method'] == 'blast_fragments':
            blast_db, db_files_to_remove = \
                build_blast_db_from_fasta_path(params['reference_seqs_fp'],
                                               output_dir=working_dir)
            self.files_to_remove += db_files_to_remove
            params['blast_db'] = blast_db
        elif params['chimera_detection_method'] == 'ChimeraSlayer':
            # copy the reference files to working dir
            # ChimeraSlayer creates an index file of the ref and
            # will crash without write permission in the ref seqs dir
            aligned_reference_seqs_fp = params['aligned_reference_seqs_fp']
            _, new_ref_filename = split(aligned_reference_seqs_fp)
            copy(aligned_reference_seqs_fp, working_dir)
            aligned_reference_seqs_fp = working_dir + "/" + new_ref_filename

            self.files_to_remove.append(aligned_reference_seqs_fp)
            params['aligned_reference_seqs_fp'] = aligned_reference_seqs_fp

            # if given, also copy the unaligned ref db
            reference_seqs_fp = params['reference_seqs_fp']
            if reference_seqs_fp:
                _, new_ref_filename = split(reference_seqs_fp)
                copy(reference_seqs_fp, working_dir)
                reference_seqs_fp = working_dir + "/" + new_ref_filename
            else:
                # otherwise create it
                reference_seqs_fp = write_degapped_fasta_to_file(
                    parse_fasta(open(aligned_reference_seqs_fp)),
                    tmp_dir=working_dir)
            # delete it afterwards
            self.files_to_remove.append(reference_seqs_fp)
            params['reference_seqs_fp'] = reference_seqs_fp

            # build blast db of reference, otherwise ChimeraSlayer will do it
            # and parallel jobs clash
            _, db_files_to_remove = \
                build_blast_db_from_fasta_path(reference_seqs_fp)
            self.files_to_remove += db_files_to_remove

            # make the index file globally
            # Reason: ChimeraSlayer first checks to see if the index file is
            # there. If not it tries to create it. This can lead to race
            # condition if several parallel jobs try to create it at the same
            # time.
            make_cidx_file(aligned_reference_seqs_fp)
            self.files_to_remove.append(aligned_reference_seqs_fp + ".cidx")
        else:
            raise ValueError("Unrecognized chimera detection method '%s'." %
                             params['chimera_detection_method'])

    def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
                          working_dir, command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate identify_chimeric_seqs.py commands which should be run."""
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_chimeric.txt']

        # Create lists to store the results.
        commands = []
        result_filepaths = []

        # Iterate over the input files.
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = \
                self._get_rename_command([fn % i for fn in out_filenames],
                                         working_dir, output_dir)
            result_filepaths += current_result_filepaths

            optional_options = ""
            if params['chimera_detection_method'] == 'blast_fragments':
                command = \
                    '%s %s -i %s -t %s -m blast_fragments -o %s -n %s -d %s -e %s -b %s %s %s' % \
                    (command_prefix,
                     self._script_name,
                     fasta_fp,
                     params['id_to_taxonomy_fp'],
                     working_dir + "/" + out_filenames[0] % i,
                     params['num_fragments'],
                     params['taxonomy_depth'],
                     params['max_e_value'],
                     params['blast_db'],
                     rename_command,
                     command_suffix)
            elif params['chimera_detection_method'] == 'ChimeraSlayer':
                optional_options = ""
                if params['min_div_ratio']:
                    optional_options += " --min_div_ratio %s" % \
                                        params['min_div_ratio']
                if params['reference_seqs_fp']:
                    optional_options += " -r %s" % params['reference_seqs_fp']

                command = \
                    '%s %s -i %s -a %s -m ChimeraSlayer -o %s %s %s %s' % \
                    (command_prefix,
                     self._script_name,
                     fasta_fp,
                     params['aligned_reference_seqs_fp'],
                     working_dir + "/" + out_filenames[0] % i,
                     optional_options,
                     rename_command,
                     command_suffix)
            else:
                raise NotImplementedError
            commands.append(command)
        return commands, result_filepaths

    def _get_poller_command(self,
                            expected_files_filepath,
                            merge_map_filepath,
                            deletion_list_filepath,
                            command_prefix='/bin/bash; ',
                            command_suffix='; exit'):
        """Generate command to initiate a poller to monitior/process completed runs
        """
        result = '%s poller.py -f %s -p %s -m %s -d %s -t %d %s' % \
            (command_prefix,
             expected_files_filepath,
             self._process_run_results_f,
             merge_map_filepath,
             deletion_list_filepath,
             self._seconds_to_sleep,
             command_suffix)
        return result, []

    def _write_merge_map_file(self, input_file_basename, job_result_filepaths,
                              params, output_dir, merge_map_filepath,
                              failures=False):
        f = open(merge_map_filepath, 'w')
        out_filepaths = [params['output_fp']]

        chims_fps = []
        logs_fps = []  # logs_fp currently not used

        for fp in job_result_filepaths:
            if fp.endswith('_chimeric.txt'):
                chims_fps.append(fp)
            else:
                log_fps.append(fp)

        for in_files, out_file in zip([chims_fps], out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()


def basic_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles

        f: file containing one set of mapping instructions per line

        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log

        If f contained the two lines above, this function would
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    infiles_lists, out_filepaths = parse_tmp_to_final_filepath_map_file(f)
    for infiles_list, out_filepath in zip(infiles_lists, out_filepaths):
        try:
            of = open(out_filepath, 'w')
        except IOError:
            raise IOError("Poller can't open final output file: %s" % out_filepath +
                          "\nLeaving individual jobs output.\n Do you have write access?")

        for fp in infiles_list:
            for line in open(fp):
                of.write('%s\n' % line.strip('\n'))
        of.close()
    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback,
    # you'll get an error right away rather than going into an infinite loop
    return True
