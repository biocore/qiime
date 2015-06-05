#!/usr/bin/env python
# File created on 25 July 201
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"


from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from qiime.parallel.pick_otus import ParallelPickOtusTrie

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Parallel pick otus using a trie"""
script_info[
    'script_description'] = """This script performs like the pick_otus.py script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel. The script uses the first p bases of each read to sort all reads into separate buckets and then each buckets is processed separately. Note that in cases of amplicon sequencing we do not expect the buckets to be even sized, but rather a few buckets make up the majority of reads. Thus, not all combination of prefix length p and number of CPUS -O make sense. Good combinations for a small desktop multicore system would be -p 5 (default) and -O 4. For larger clusters, we suggest -p 10 and -O 20. Increasing -p to a value much larger than 10 will lead to lots of temporary files and many small jobs, so likely will not speed up the OTU picking. On the other hand, the max speed-up is bounded by the size of the largest buckets, so adding more cores will not always increase efficiency."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example""",
     """Pick OTUs by building a trie out of $PWD/inseqs.fasta and write the output to the $PWD/trie_otus/ directory. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).""",
     """%prog -i $PWD/seqs.fna -o $PWD/trie_otus/"""))
script_info['script_usage'].append(
    ("""Example""",
     """Pick OTUs by building a trie out of $PWD/inseqs.fasta and write the output to the $PWD/trie_otus/ directory. Split the input according to the first 10 bases of each read and process each set independently.""",
     """%prog -i $PWD/seqs.fna -o $PWD/trie_otus/ -p 10"""))

script_info[
    'output_description'] = """The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file."""

script_info['required_options'] = [
    make_option('-i', '--input_fasta_fp', action='store',
                type='existing_filepath', help='full path to ' +
                'input_fasta_fp'),
    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', help='path to store output files'),
]

script_info['optional_options'] = [
    make_option('-p', '--prefix_length', action='store',
                type='int', help='prefix length used to split the input. '
                +
                'Must be smaller than the shortest seq in input! [default: %default]',
                default=5),
    options_lookup['jobs_to_start'],
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['poll_directly'],
    options_lookup['cluster_jobs_fp'],
    options_lookup['suppress_polling'],
    options_lookup['job_prefix'],
    options_lookup['seconds_to_sleep']
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # create dict of command-line options
    params = eval(str(opts))

    parallel_runner = ParallelPickOtusTrie(
        cluster_jobs_fp=opts.cluster_jobs_fp,
        jobs_to_start=opts.jobs_to_start,
        retain_temp_files=opts.retain_temp_files,
        suppress_polling=opts.suppress_polling,
        seconds_to_sleep=opts.seconds_to_sleep)
    parallel_runner(opts.input_fasta_fp,
                    opts.output_dir,
                    params,
                    job_prefix=opts.job_prefix,
                    poll_directly=opts.poll_directly,
                    suppress_submit_jobs=opts.suppress_submit_jobs)

if __name__ == "__main__":
    main()
