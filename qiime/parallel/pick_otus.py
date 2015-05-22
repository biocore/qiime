#!/usr/bin/env python
# File created on 07 Jul 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jens Reeder", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from math import ceil
from os.path import basename, join
from os import mkdir
from re import compile

from bfillings.formatdb import build_blast_db_from_fasta_path
from bfillings.sortmerna_v2 import build_database_sortmerna

from skbio.parse.sequences import parse_fasta

from qiime.parallel.util import ParallelWrapper, BufferedWriter
from qiime.parallel.poller import basic_process_run_results_f


class ParallelPickOtus(ParallelWrapper):
    _script_name = "pick_otus.py"
    _job_prefix = 'POTU'
    _input_splitter = ParallelWrapper._split_fasta
    _process_run_results_f =\
        'qiime.parallel.pick_otus.parallel_pick_otus_process_run_results_f'

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

    def _write_merge_map_file(self,
                              input_file_basename,
                              job_result_filepaths,
                              params,
                              output_dir,
                              merge_map_filepath):
        """
        """
        f = open(merge_map_filepath, 'w')

        otus_fps = []
        log_fps = []
        failures_fps = []

        # determine if the otu picking method generates
        # failures files
        failures = False
        for fp in job_result_filepaths:
            if fp.endswith('_failures.txt'):
                failures = True
                break

        if not failures:
            out_filepaths = [
                '%s/%s_otus.txt' % (output_dir, input_file_basename),
                '%s/%s_otus.log' % (output_dir, input_file_basename)]
            in_filepaths = [otus_fps, log_fps]
        else:
            out_filepaths = [
                '%s/%s_otus.txt' % (output_dir, input_file_basename),
                '%s/%s_otus.log' % (output_dir, input_file_basename),
                '%s/%s_failures.txt' % (output_dir, input_file_basename)]
            in_filepaths = [otus_fps, log_fps, failures_fps]

        for fp in job_result_filepaths:
            if fp.endswith('_otus.txt'):
                otus_fps.append(fp)
            elif fp.endswith('_otus.log'):
                log_fps.append(fp)
            elif fp.endswith('_failures.txt'):
                failures_fps.append(fp)
            else:
                pass

        for in_files, out_file in\
                zip(in_filepaths, out_filepaths):
            f.write('\t'.join(in_files + [out_file]))
            f.write('\n')
        f.close()


class ParallelPickOtusSortMeRNA(ParallelPickOtus):
    """ Run pick_otus.py in parallel for SortMeRNA
    """

    def _precommand_initiation(
            self, input_fp, output_dir, working_dir, params):
        if not params['sortmerna_db']:
            # Build the sortmerna database from the reference_seqs_fp
            sortmerna_db, db_files_to_remove = \
                build_database_sortmerna(params['refseqs_fp'],
                                         max_pos=params['sortmerna_max_pos'],
                                         output_dir=working_dir)
            self.files_to_remove += db_files_to_remove
            params['sortmerna_db'] = sortmerna_db

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be submitted to cluster
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_otus.log',
                         job_prefix + '.%d_otus.txt',
                         job_prefix + '.%d_failures.txt']

        # Create lists to store the results
        commands = []
        result_filepaths = []

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # temporary working directory for sortmerna per job                                                                                                                                                                                                                         
            working_dir_t = join(working_dir, "%d" % i)
            mkdir(working_dir_t)

            # Each run ends with moving the output file from the tmp dir to                                                                                                                                                                                                              
            # the output_dir. Build the command to perform the move here.                                                                                                                                                                                                                
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % i for fn in out_filenames], working_dir_t, output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -m sortmerna -i %s -r %s --sortmerna_db %s -o %s --sortmerna_e_value %s -s %s --threads %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 params['refseqs_fp'],
                 params['sortmerna_db'],
                 working_dir_t,
                 params['sortmerna_e_value'],
                 params['similarity'],
                 params['threads'],
                 rename_command,
                 command_suffix)

            commands.append(command)

        return commands, result_filepaths


class ParallelPickOtusUclustRef(ParallelPickOtus):

    def _identify_files_to_remove(self, job_result_filepaths, params):
        """ Select the files to remove: by default remove all files
        """
        if params['save_uc_files']:
            # keep any .uc files that get created
            result =\
                [fp for fp in job_result_filepaths if not fp.endswith('.uc')]
        else:
            result = [job_result_filepaths]

        return result

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be run
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_otus.log',
                         job_prefix + '.%d_otus.txt',
                         job_prefix + '.%s_failures.txt']

        # Create lists to store the results
        commands = []
        result_filepaths = []

        if params['enable_rev_strand_match']:
            enable_rev_strand_match_str = '-z'
        else:
            enable_rev_strand_match_str = ''

        if params['optimal_uclust']:
            optimal_uclust_str = '-A'
        else:
            optimal_uclust_str = ''

        if params['exact_uclust']:
            exact_uclust_str = '-E'
        else:
            exact_uclust_str = ''

        if params['stable_sort']:
            stable_sort_str = ''
        else:
            stable_sort_str = '--suppress_uclust_stable_sort'

        if params['save_uc_files']:
            save_uc_files_str = ''
            out_filenames += [job_prefix + '.%d_clusters.uc']
        else:
            save_uc_files_str = '-d'

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % i for fn in out_filenames],
                working_dir,
                output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -i %s -r %s -m uclust_ref --suppress_new_clusters -o %s -s %s %s %s %s --max_accepts %s --max_rejects %s --stepwords %d --w %d %s %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 params['refseqs_fp'],
                 working_dir,
                 params['similarity'],
                 enable_rev_strand_match_str,
                 optimal_uclust_str,
                 exact_uclust_str,
                 params['max_accepts'],
                 params['max_rejects'],
                 params['stepwords'],
                 params['word_length'],
                 stable_sort_str,
                 save_uc_files_str,
                 rename_command,
                 command_suffix)

            commands.append(command)

        return commands, result_filepaths


class ParallelPickOtusUsearch61Ref(ParallelPickOtus):

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be run
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_otus.log',
                         job_prefix + '.%d_otus.txt',
                         job_prefix + '.%s_failures.txt']

        # Create lists to store the results
        commands = []
        result_filepaths = []

        # Generate the parameters to pass to pick_otus.py. This must exclude
        # parameters that get passed only to the parallel version
        # (e.g jobs_to_start) and values that get overwritten (e.g.,
        # input_fasta_fp)
        param_fields = []
        ignored_params = set(["input_fasta_fp", "output_dir", "jobs_to_start",
                              "retain_temp_files", "suppress_submit_jobs", "poll_directly",
                              "cluster_jobs_fp", "suppress_polling", "job_prefix",
                              "seconds_to_sleep"])
        for name, value in params.items():
            if name in ignored_params or value == False:
                pass
            elif value == 'True':
                param_fields.append('--%s' % name)
            else:
                param_fields.append('--%s %s' % (name, value))
        params_str = ' '.join(param_fields)

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            iteration_working_dir = join(working_dir, str(i))
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % i for fn in out_filenames],
                iteration_working_dir,
                output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s %s -i %s -m usearch61_ref --suppress_new_clusters -o %s %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 iteration_working_dir,
                 params_str,
                 rename_command,
                 command_suffix)
            commands.append(command)

        return commands, result_filepaths


class ParallelPickOtusBlast(ParallelPickOtus):

    def _precommand_initiation(
            self, input_fp, output_dir, working_dir, params):
        if not params['blast_db']:
            # Build the blast database from the reference_seqs_fp -- all procs
            # will then access one db rather than create one per proc
            blast_db, db_files_to_remove = \
                build_blast_db_from_fasta_path(params['refseqs_fp'])
            self.files_to_remove += db_files_to_remove
            params['blast_db'] = blast_db

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be submitted to cluster
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = [job_prefix + '.%d_otus.log',
                         job_prefix + '.%d_otus.txt']

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
                '%s %s -i %s -b %s -m blast -o %s -e %s -s %s --min_aligned_percent %s %s %s' %\
                (command_prefix,
                 self._script_name,
                 fasta_fp,
                 params['blast_db'],
                 working_dir,
                 params['max_e_value'],
                 params['similarity'],
                 params['min_aligned_percent'],
                 rename_command,
                 command_suffix)

            commands.append(command)

        return commands, result_filepaths


def parallel_pick_otus_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles

        f: file containing one set of mapping instructions per line

        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log
         f1_failures.txt f2_failures.txt f3_failures.txt f_failires.txt

        If f contained the two lines above, this function would
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log
    """
    lines = list(f)
    # handle catting of log files and failure files
    basic_process_run_results_f([lines[1]])
    try:
        basic_process_run_results_f([lines[2]])
        basic_process_run_results_f([lines[3]])
    except IndexError:
        # no failures files or blast6 were generated (BLAST
        # doesn't create these)
        pass
    # handle merging of otu maps
    fields = lines[0].strip().split()
    infiles_list = fields[:-1]
    out_filepath = fields[-1]
    try:
        of = open(out_filepath, 'w')
    except IOError:
        raise IOError("Poller can't open final output file: %s" % out_filepath +
                      "\nLeaving individual jobs output.\n Do you have write access?")

    unique_otu_map = {}
    for fp in infiles_list:
        for line in open(fp):
            fields = line.strip().split()
            try:
                # current otu_id already exists, so append this
                # set of seq_ids
                unique_otu_map[fields[0]] += fields[1:]
            except KeyError:
                # current otu_id has not been seen yet, so
                # create it with the current set of otus
                unique_otu_map[fields[0]] = fields[1:]

    for otu_id, seq_ids in unique_otu_map.items():
        of.write('\t'.join([otu_id] + seq_ids))
        of.write('\n')
    of.close()

    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback,
    # you'll get an error right away rather than going into an infinite loop
    return True


class ParallelPickOtusTrie(ParallelPickOtus):

    """Picking Otus using a trie the parallel way

    We parallelize the Trie OTU picker using this scheme:
    1. use the exact prefix filter with a short wordlength (say 5)
       to sort all reads into buckets according to their first 5 nucs.
    2. Run trie otupicker an each bucket separately, distribute over cluster
    3. Combine mappings of 2, Since each bucket is independent fro the rest,
       a simple cat with incrementing OTU ids should do it.
    """
    _process_run_results_f =\
        'qiime.parallel.pick_otus.parallel_pick_otus_trie_process_run_results_f'

    def _call_initialization(self,
                             input_fp,
                             output_dir,
                             params,
                             job_prefix,
                             poll_directly,
                             suppress_submit_jobs):
        """ Called as the first step in __call__.
        """
        self.prefix_counts = {}

    def _split_along_prefix(self,
                            input_fp,
                            params,
                            jobs_to_start,
                            job_prefix,
                            output_dir):
        """ Split input sequences into sets with identical prefix"""
        out_files = []
        buffered_handles = {}
        prefix_length = params['prefix_length'] or 1
        for seq_id, seq in parse_fasta(open(input_fp)):

            if(len(seq) < prefix_length):
                raise ValueError("Prefix length must be equal or longer than sequence.\n"
                                 + " Found seq %s with length %d" % (seq_id, len(seq)))
            prefix = seq[:prefix_length]

            if (prefix not in buffered_handles):
                # never seen this prefix before
                out_fp = "%s/%s%s" % (output_dir, job_prefix, prefix)
                buffered_handles[prefix] = BufferedWriter(out_fp)
                out_files.append(out_fp)
                self.prefix_counts[prefix] = 0

            self.prefix_counts[prefix] += 1
            buffered_handles[prefix].write('>%s\n%s\n' % (seq_id, seq))

        # make sure all buffers are closed and flushed
        for buf_fh in buffered_handles.itervalues():
            buf_fh.close()

        remove_files = True
        return out_files, remove_files

    _input_splitter = _split_along_prefix

    def _get_job_commands(self,
                          fasta_fps,
                          output_dir,
                          params,
                          job_prefix,
                          working_dir,
                          command_prefix='/bin/bash; ',
                          command_suffix='; exit'):
        """Generate pick_otus commands which should be run
        """
        # Create basenames for each of the output files. These will be filled
        # in to create the full list of files created by all of the runs.
        out_filenames = ['%s_otus.log',
                         '%s_otus.txt']

        # Create lists to store the results
        commands = {}
        result_filepaths = []

        # Iterate over the input files
        for i, fasta_fp in enumerate(fasta_fps):
            # Each run ends with moving the output file from the tmp dir to
            # the output_dir. Build the command to perform the move here.
            rename_command, current_result_filepaths = self._get_rename_command(
                [fn % basename(fasta_fp) for fn in out_filenames],
                working_dir,
                output_dir)
            result_filepaths += current_result_filepaths

            command = \
                '%s -i %s -m trie -o %s %s' %\
                (self._script_name,
                 fasta_fp,
                 working_dir,
                 rename_command)

            re = compile("POTU_\w+_(\w+)")
            prefix = (re.search(fasta_fp)).group(1)
            commands[prefix] = command

        commands = self._merge_to_n_commands(commands,
                                             params['jobs_to_start'],
                                             command_prefix=command_prefix,
                                             command_suffix=command_suffix)
        return commands, result_filepaths

    def _merge_to_n_commands(self,
                             commands,
                             n,
                             delimiter=' ; ',
                             command_prefix=None,
                             command_suffix=None):
        """ merge a list of commands into n commands

        Uses the size of each jobs to estimate an even distribution
        of commands ofver the n jobs.

        commands: dict of commands keyed by prefix

        n: number of jobs
        """
        if n < 1:
            raise ValueError("number of commands (n) must be an integer >= 1")

        # Distribute jobs according to their load
        grouped_prefixe, levels = greedy_partition(self.prefix_counts, n)

        # TODO: remove after profiling or maybe move to log file
#       print levels
#       print ("Maximal theoretical speed up: %f" % (sum(levels)/max(levels)))

        result = []
        # begin iterating through the commands
        current_cmds = []

        for bucket in grouped_prefixe:
            for prefix in bucket:
                command = commands[prefix]

                subcommands = [c.strip() for c in command.split(';')]
                current_cmds.append(delimiter.join(subcommands))

            result.append(delimiter.join(current_cmds))
            current_cmds = []

        for i, r in enumerate(result):
            r = '%s %s %s' % (command_prefix, r, command_suffix)
            result[i] = r.strip()

        return result


def parallel_pick_otus_trie_process_run_results_f(f):
    """ Copy each list of infiles to each outfile and delete infiles

        f: file containing one set of mapping instructions per line

        example f:
         f1.txt f2.txt f3.txt f_combined.txt
         f1.log f2.log f3.log f_combined.log
         f1_failures.txt f2_failures.txt f3_failures.txt f_failires.txt

         If f contained the two lines above, this function would
         concatenate f1.txt, f2.txt, and f3.txt into f_combined.txt
         and f1.log, f2.log, and f3.log into f_combined.log

         Note: this is mostly copied from
               parallel_pick_otus_process_run_results_f
               The only difference is the way how otu_maps are combined.
               Since each otu_map for the parallel trie otu pickers yields disjoint
               sets, we have to make the otu_ids disjoint as well.
    """
    lines = list(f)
    # handle catting of log files and failure files
    basic_process_run_results_f([lines[1]])
    try:
        basic_process_run_results_f([lines[2]])
        basic_process_run_results_f([lines[3]])
    except IndexError:
        # no failures files or blast6 were generated (BLAST
        # doesn't create these)
        pass
    # handle merging of otu maps
    fields = lines[0].strip().split()
    infiles_list = fields[:-1]
    out_filepath = fields[-1]
    try:
        of = open(out_filepath, 'w')
    except IOError:
        raise IOError("Poller can't open final output file: %s" % out_filepath +
                      "\nLeaving individual jobs output.\n Do you have write access?")

    unique_otu_map = {}
    otu_id = 0
    for fp in infiles_list:
        for line in open(fp):
            fields = line.strip().split()
            of.write('\t'.join(["%d" % otu_id] + fields[1:]))
            of.write('\n')
            otu_id += 1
    of.close()

    # It is a good idea to have your clean_up_callback return True.
    # That way, if you get mixed up and pass it as check_run_complete_callback,
    # you'll get an error right away rather than going into an infinite loop
    return True


def greedy_partition(counts, n):
    """Distribute k counts evenly across n buckets,

    counts: dict of key, counts pairs
    n: number of buckets that the counts should be distributed over
    """

    buckets = [[] for i in range(n)]
    fill_levels = [0 for i in range(n)]

    for key in sorted(counts, reverse=True,
                      key=lambda c: counts[c]):
        smallest = fill_levels.index(min(fill_levels))
        buckets[smallest].append(key)
        fill_levels[smallest] += counts[key]

    return buckets, fill_levels
