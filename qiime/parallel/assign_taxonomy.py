#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from bfillings.formatdb import build_blast_db_from_fasta_path

from qiime.parallel.util import ParallelWrapper


class ParallelTaxonomyAssigner(ParallelWrapper):
    _script_name = 'assign_taxonomy.py'
    _input_splitter = ParallelWrapper._split_fasta

    def _build_job_commands(self, tax_specific_param_str, fasta_fps,
                            output_dir, params, job_prefix, working_dir,
                            command_prefix='/bin/bash; ',
                            command_suffix='; exit'):
        """Generate assign_taxonomy.py commands which should be run."""
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_tax_assignments.log',
                         job_prefix + '.%d_tax_assignments.txt']

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

            command = '%s %s %s -o %s -i %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 tax_specific_param_str,
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
                              merge_map_filepath,
                              failures=False):
        """
        """
        f = open(merge_map_filepath, 'w')

        out_filepaths = [
            '%s/%s_tax_assignments.txt' % (output_dir, input_file_basename),
            '%s/%s_tax_assignments.log' % (output_dir, input_file_basename)]

        assignment_fps = []
        log_fps = []

        for fp in job_result_filepaths:
            if fp.endswith('_tax_assignments.txt'):
                assignment_fps.append(fp)
            else:
                log_fps.append(fp)

        for in_files, out_file in\
                zip([assignment_fps, log_fps], out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()


class ParallelRdpTaxonomyAssigner(ParallelTaxonomyAssigner):
    _job_prefix = 'RDP'

    def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
                          working_dir, command_prefix=None,
                          command_suffix='; exit'):
        command_prefix = command_prefix or \
            '/bin/bash; export RDP_JAR_PATH=%s; ' % params['rdp_classifier_fp']

        rdp_params = '-m rdp -c %1.2f --rdp_max_memory %d ' % (
            params['confidence'], params['rdp_max_memory'])
        if params['id_to_taxonomy_fp'] and params['reference_seqs_fp']:
            rdp_params += '-t %s -r %s' % (params['id_to_taxonomy_fp'],
                                           params['reference_seqs_fp'])

        return self._build_job_commands(rdp_params, fasta_fps, output_dir,
                                        params, job_prefix, working_dir, command_prefix,
                                        command_suffix)


class ParallelUclustConsensusTaxonomyAssigner(ParallelTaxonomyAssigner):
    _job_prefix = 'UCTA'

    def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
                          working_dir, command_prefix=None,
                          command_suffix='; exit'):
        command_prefix = command_prefix or ''

        uclust_params = ' '.join(
            ['-m uclust',
             '--min_consensus_fraction %f' % params['min_consensus_fraction'],
             '--similarity %f' % params['similarity'],
             '--uclust_max_accepts %d' % params['uclust_max_accepts'],
             '-t %s' % params['id_to_taxonomy_fp'],
             '-r %s' % params['reference_seqs_fp']])

        return self._build_job_commands(uclust_params,
                                        fasta_fps,
                                        output_dir,
                                        params,
                                        job_prefix,
                                        working_dir,
                                        command_prefix,
                                        command_suffix)


class ParallelBlastTaxonomyAssigner(ParallelTaxonomyAssigner):
    _job_prefix = 'BTA'

    def _precommand_initiation(
            self, input_fp, output_dir, working_dir, params):
        if not params['blast_db']:
            # Build the blast database from the reference_seqs_fp -- all procs
            # will then access one db rather than create one per proc.
            blast_db, db_files_to_remove = build_blast_db_from_fasta_path(
                params['reference_seqs_fp'], output_dir=working_dir)
            self.files_to_remove += db_files_to_remove
            params['blast_db'] = blast_db

    def _get_job_commands(self, fasta_fps, output_dir, params, job_prefix,
                          working_dir, command_prefix=None,
                          command_suffix='; exit'):
        command_prefix = command_prefix or \
            '/bin/bash; cd %s; export BLASTMAT=%s;' % (working_dir,
                                                       params['blastmat_dir'])

        blast_params = '-m blast -e %s -b %s -t %s ' % (
            params['e_value'], params['blast_db'],
            params['id_to_taxonomy_fp'])

        return self._build_job_commands(blast_params, fasta_fps, output_dir,
                                        params, job_prefix, working_dir, command_prefix,
                                        command_suffix)
