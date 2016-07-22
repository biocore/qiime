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

from os.path import join

from qiime.make_otu_table import make_otu_table
from qiime.parse import parse_observation_metadata
from qiime.parallel.pick_otus import ParallelPickOtus
from qiime.util import write_biom_table


class ParallelDatabaseMapper(ParallelPickOtus):
    _script_name = 'map_reads_to_reference.py'
    _job_prefix = 'RMAP'

    def _call_cleanup(self,
                      input_fp,
                      output_dir,
                      params,
                      job_prefix,
                      poll_directly,
                      suppress_submit_jobs):
        """ Called as the last step in __call__.
        """
        if poll_directly:
            if params['observation_metadata_fp'] is not None:
                observation_metadata = \
                    parse_observation_metadata(
                        open(params['observation_metadata_fp'], 'U'))
            else:
                observation_metadata = None
            biom_fp = join(output_dir, 'observation_table.biom')
            biom_table = make_otu_table(
                open(join(output_dir, 'observation_map.txt'), 'U'),
                observation_metadata)
            write_biom_table(biom_table, biom_fp)
        else:
            # can't construct the final biom file if not polling
            # directly as the final observation map won't have been created yet
            pass


class ParallelDatabaseMapperUsearch(ParallelDatabaseMapper):

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        out_filenames = ['observation_map.txt',
                         'out.uc',
                         'out.bl6',
                         'observation_table.biom']

        # Create lists to store the results
        commands = []
        result_filepaths = []

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            run_output_dir = join(working_dir, str(i))
            tmp_output_dir = join(working_dir, str(i), 'tmp')
            rename_command, current_result_filepaths = self._get_rename_command(
                out_filenames,
                tmp_output_dir,
                run_output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -i %s -r %s -m usearch -o %s --min_percent_id %s --max_accepts %d --max_rejects %d --queryalnfract %f --targetalnfract %f --evalue %e %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 params['refseqs_fp'],
                 tmp_output_dir,
                 params['min_percent_id'],
                 params['max_accepts'],
                 params['max_rejects'],
                 params['queryalnfract'],
                 params['targetalnfract'],
                 params['evalue'],
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
        """
        """
        f = open(merge_map_filepath, 'w')

        observation_fps = []
        uc_fps = []
        blast6_fps = []

        out_filepaths = [
            '%s/observation_map.txt' % output_dir,
            '%s/out.uc' % output_dir,
            '%s/out.bl6' % output_dir]
        in_filepaths = [observation_fps, uc_fps, blast6_fps]

        for fp in job_result_filepaths:
            if fp.endswith('observation_map.txt'):
                observation_fps.append(fp)
            if fp.endswith('.uc'):
                uc_fps.append(fp)
            elif fp.endswith('.bl6'):
                blast6_fps.append(fp)
            else:
                pass

        for in_files, out_file in\
                zip(in_filepaths, out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()


class ParallelDatabaseMapperBlat(ParallelDatabaseMapper):

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        out_filenames = ['observation_map.txt',
                         'out.bl9',
                         'observation_table.log',
                         'observation_table.biom']

        # Create lists to store the results
        commands = []
        result_filepaths = []
        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            run_output_dir = join(working_dir, str(i))
            tmp_output_dir = join(working_dir, str(i), 'tmp')
            rename_command, current_result_filepaths = self._get_rename_command(
                out_filenames,
                tmp_output_dir,
                run_output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -i %s -r %s -m blat -o %s --min_percent_id %s --evalue %e %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 params['refseqs_fp'],
                 tmp_output_dir,
                 params['min_percent_id'],
                 params['evalue'],
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
        """
        """
        f = open(merge_map_filepath, 'w')

        observation_fps = []
        log_fps = []
        blast9_fps = []

        out_filepaths = [
            '%s/observation_map.txt' % output_dir,
            '%s/observation_table.log' % output_dir,
            '%s/out.bl9' % output_dir]
        in_filepaths = [observation_fps, log_fps, blast9_fps]

        for fp in job_result_filepaths:
            if fp.endswith('observation_map.txt'):
                observation_fps.append(fp)
            if fp.endswith('.log'):
                log_fps.append(fp)
            elif fp.endswith('.bl9'):
                blast9_fps.append(fp)
            else:
                pass

        for in_files, out_file in\
                zip(in_filepaths, out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()


class ParallelDatabaseMapperBwaShort(ParallelDatabaseMapper):

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        out_filenames = ['observation_map.txt',
                         'bwa_raw_out.sam',
                         'bwa_raw_out.sai',
                         'observation_table.log',
                         'observation_table.biom']

        # Create lists to store the results
        commands = []
        result_filepaths = []
        if params['max_diff'] is not None:
            max_diff_str = "--max_diff %s" % params['max_diff']
        else:
            max_diff_str = ""

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            run_output_dir = join(working_dir, str(i))
            tmp_output_dir = join(working_dir, str(i), 'tmp')
            rename_command, current_result_filepaths = self._get_rename_command(
                out_filenames,
                tmp_output_dir,
                run_output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -i %s -r %s -m bwa-short -o %s %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 params['refseqs_fp'],
                 tmp_output_dir,
                 max_diff_str,
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
        """
        """
        f = open(merge_map_filepath, 'w')

        observation_fps = []
        log_fps = []
        sam_fps = []

        out_filepaths = [
            '%s/observation_map.txt' % output_dir,
            '%s/observation_table.log' % output_dir,
            '%s/bwa_raw_out.sam' % output_dir]
        in_filepaths = [observation_fps, log_fps, sam_fps]

        for fp in job_result_filepaths:
            if fp.endswith('observation_map.txt'):
                observation_fps.append(fp)
            if fp.endswith('.log'):
                log_fps.append(fp)
            elif fp.endswith('.sam'):
                sam_fps.append(fp)
            else:
                pass

        for in_files, out_file in\
                zip(in_filepaths, out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()
