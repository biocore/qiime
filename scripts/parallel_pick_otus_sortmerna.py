#!/usr/bin/env python

from __future__ import division
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jenya Kopylov"
__email__ = "jenya.kopylov@gmail.com"

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from qiime.parallel.pick_otus import ParallelPickOtusSortMeRNA

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Parallel pick otus using SortMeRNA"
script_info['script_description'] = (
    "This script works like the pick_otus.py script, but is intended to make "
    "use of multicore/multiprocessor environments to perform analyses in "
    "parallel.")
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Example",
     "Pick OTUs by searching $PWD/inseqs.fasta against $PWD/refseqs.fasta "
     "with reference-based sortmerna and write the output to the "
     "$PWD/smr_otus/ directory. This is a closed-reference OTU picking "
     "process. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented "
     "here as $PWD, but will generally look something like "
     "/home/ubuntu/my_analysis/).",
     "%prog -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/smr_otus/"))
script_info['output_description'] = (
    "The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). "
    "The .txt file is composed of tab-delimited lines, where the first field "
    "on each line corresponds to an OTU identifier which is the reference "
    "sequence identifier, and the remaining fields correspond to sequence "
    "identifiers assigned to that OTU. The resulting .log file contains a "
    "list of parameters passed to this script along with the output location "
    "of the resulting .txt file.")

script_info['required_options'] = [
    make_option('-i', '--input_fasta_fp', action='store',
                type='existing_filepath', help='Path to input fasta file.'),

    make_option('-o', '--output_dir', action='store',
                type='new_dirpath', help='Directory where output should be written.'),

    make_option('-r', '--refseqs_fp', action='store',
                type='existing_filepath', help='Path to reference fasta file.')
]

script_info['optional_options'] = [

    make_option('-s', '--similarity', action='store',
                type='float', help='Sequence similarity ' +
                'threshold [default: %default]', default=0.97),

    make_option('--sortmerna_db', type='string',
                help='Pre-existing database to search against when using '
                     '-m sortmerna [default: %default]'),

    make_option('--sortmerna_e_value', type='float', default=1,
                help='Maximum E-value when clustering [default = %default]'),

    make_option('--sortmerna_coverage', type='float', default=0.97,
                help='Mininum percent query coverage (of an alignment) '
                     'to consider a hit, expressed as a fraction between 0 '
                     'and 1 [default: %default]'),

    make_option('--sortmerna_tabular', default=False, action='store_true',
                help='Output alignments in the Blast-like tabular format '
                     'with two additional columns including the CIGAR '
                     'string and the percent query coverage '
                     '[default: %default]'),

    make_option('--sortmerna_best_N_alignments', type='int', default=1,
                help='Must be set together with --sortmerna_tabular. '
                     'This option specifies how many alignments per read '
                     'will be written [default: %default]'),

    make_option('--sortmerna_max_pos', type='int', default=10000,
                help='The maximum number of positions per seed to store '
                     ' in the indexed database [default: %default]'),

    make_option('--threads', default=1, help=
                "Specify the number of threads to use per job. "
                "Use --jobs_to_start to specify the number of jobs."
                "[default: %default]"),

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

    parallel_runner = ParallelPickOtusSortMeRNA(
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
                    suppress_submit_jobs=False)


if __name__ == "__main__":
    main()
