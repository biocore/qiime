#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from bfillings.formatdb import build_blast_db_from_fasta_path

from qiime.align_seqs import compute_min_alignment_length
from qiime.parallel.util import ParallelWrapper
from qiime.util import get_qiime_temp_dir


class ParallelAlignSeqsPyNast(ParallelWrapper):
    _script_name = "align_seqs.py"
    _job_prefix = 'ALIGN'
    _input_splitter = ParallelWrapper._split_fasta

    def _precommand_initiation(
            self, input_fp, output_dir, working_dir, params):
        if params['min_length'] < 0:
            params['min_length'] = compute_min_alignment_length(
                open(input_fp, 'U'))

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate PyNAST commands which should be submitted to cluster
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_aligned.fasta',
                         job_prefix + '.%d_failures.fasta',
                         job_prefix + '.%d_log.txt']

        # Initialize the command_prefix and command_suffix
        command_prefix = command_prefix or '/bin/bash; '
        command_suffix = command_suffix or '; exit'

        # Create lists to store the results
        commands = []
        result_filepaths = []

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % i for fn in out_filenames], working_dir, output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -p %1.2f -e %d -m pynast -t %s -a %s -o %s -i %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 params['min_percent_id'],
                 params['min_length'],
                 params['template_fp'],
                 params['pairwise_alignment_method'],
                 working_dir,
                 fasta_fp,
                 rename_command,
                 command_suffix)

            commands.append(command)

        return commands, result_filepaths

    def _write_merge_map_file(self,
                              input_file_basename,
                              job_result_filepaths,
                              params,
                              output_dir,
                              merge_map_filepath):

        f = open(merge_map_filepath, 'w')

        out_filepaths = [
            '%s/%s_aligned.fasta' % (output_dir, input_file_basename),
            '%s/%s_failures.fasta' % (output_dir,
                                      input_file_basename),
            '%s/%s_log.txt' % (output_dir, input_file_basename)]

        aligned_fps = []
        failures_fps = []
        log_fps = []

        for fp in job_result_filepaths:
            if fp.endswith('_aligned.fasta'):
                aligned_fps.append(fp)
            elif fp.endswith('_failures.fasta'):
                failures_fps.append(fp)
            else:
                log_fps.append(fp)

        for in_files, out_file in\
                zip([aligned_fps, failures_fps, log_fps], out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()
