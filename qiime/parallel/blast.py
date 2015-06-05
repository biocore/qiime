#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import split, splitext

from bfillings.formatdb import build_blast_db_from_fasta_path

from qiime.util import load_qiime_config, get_options_lookup
from qiime.parallel.util import ParallelWrapper


class ParallelBlaster(ParallelWrapper):
    _script_name = load_qiime_config()['blastall_fp']
    _input_splitter = ParallelWrapper._split_fasta
    _job_prefix = 'BLAST'

    def _precommand_initiation(
            self, input_fp, output_dir, working_dir, params):
        if params['refseqs_path']:
            # Build the blast database from the refseqs_path -- all procs
            # will then access one db rather than create one per proc.
            blast_db, db_files_to_remove = \
                build_blast_db_from_fasta_path(params['refseqs_path'])
            self.files_to_remove += db_files_to_remove
            params['blast_db'] = blast_db

    def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
                          working_dir, command_prefix=None,
                          command_suffix='; exit'):
        """Generate blastall commands which should be run."""
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_blast_out.txt']

        command_prefix = command_prefix or \
            '/bin/bash; export BLASTMAT=%s;' % params['blastmat_dir']

        if not params['disable_low_complexity_filter']:
            complexity_filter_str = 'T'
        else:
            complexity_filter_str = 'F'

        # Create lists to store the results.
        commands = []
        result_filepaths = []

        # Iterate over the input files.
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            # rename_command, current_result_filepaths = \
            #        self._get_rename_command([fn % i for fn in out_filenames],
            #                                 working_dir, output_dir)
            #result_filepaths += current_result_filepaths

            # TODO should this be put in self._get_rename_command()?
            infile_basename = splitext(split(fasta_fp)[1])[0]
            working_outfile_path = '%s/%s_blast_out.txt' %\
                (working_dir, infile_basename)
            outfile_path = '%s/%s_blast_out.txt' % (output_dir,
                                                    infile_basename)
            rename_command = '; mv %s %s' % (working_outfile_path,
                                             outfile_path)
            result_filepaths.append(outfile_path)

            command = '%s %s -p blastn -m 9 -e %s -F %s -W %s -b %s -i %s -d %s > %s %s %s' % \
                (command_prefix,
                 self._script_name,
                 params['e_value'],
                 complexity_filter_str,
                 params['word_size'],
                 params['num_hits'],
                 fasta_fp,
                 params['blast_db'],
                 working_outfile_path,
                 rename_command,
                 command_suffix)
            commands.append(command)
        return commands, result_filepaths

    def _write_merge_map_file(self,
                              input_file_basename,
                              job_result_filepaths,
                              params,
                              output_dir,
                              merge_map_filepath,
                              failures=False):
        """
        """
        f = open(merge_map_filepath, 'w')
        out_filepath = '%s/%s_blast_out.txt' % (output_dir,
                                                input_file_basename)
        f.write('\t'.join(job_result_filepaths + [out_filepath]))
        f.write('\n')
        f.close()
