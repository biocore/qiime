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

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from qiime.parallel.map_reads_to_reference import \
    (ParallelDatabaseMapperUsearch, ParallelDatabaseMapperBlat,
     ParallelDatabaseMapperBwaShort)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info[
    'script_usage'] = [("", "", "%prog -i $PWD/query_nt.fasta -r $PWD/refseqs_pr.fasta -o $PWD/usearch_mapped"),
                       ("", "",
                        "%prog -i $PWD/query_nt.fasta -r $PWD/refseqs_pr.fasta -o $PWD/blat_mapped -m blat"),
                       ("", "", "%prog -i $PWD/query_nt.fasta -r $PWD/refseqs_nt.fasta -o $PWD/bwa-short_mapped -m bwa-short")]
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_seqs_filepath', type='existing_filepath',
                help='Path to input sequences file'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='Directory to store results'),
    make_option('-r', '--refseqs_fp', type='existing_filepath',
                help='Path to reference sequences'),
]
script_info['optional_options'] = [
    make_option('-t', '--observation_metadata_fp', type='existing_filepath',
                help=(
                    'Path to observation metadata (e.g., taxonomy, EC, etc) [default: %default]')),

    make_option('-m', '--assignment_method', type='choice',
                choices=['usearch', 'blat', 'bwa-short'], default = "usearch",
                help=('Method for picking OTUs.  Valid choices are: ' +
                      ' '.join(['usearch', 'blat', 'bwa-short']) +
                      '. [default: %default]')),

    make_option('-e', '--evalue', type='float', default=1e-10,
                help=('Max e-value to consider a match [default: %default]')),

    make_option('-s', '--min_percent_id', type='float', default=0.75,
                help=(
                    'Min percent id to consider a match, expressed as a '
                    'fraction between 0 and 1 [default: %default]')),

    make_option('--max_diff', type='float', default=None,
                help=('maxDiff to consider a match (applicable for -m bwa) -- '
                      'see the aln section of "man bwa" for details '
                      '[default (defined by bwa): 0.04]')),

    make_option('--queryalnfract', type='float', default=0.35,
                help=(
                    'Min percent of the query seq that must match to consider '
                    'a match, expressed as a fraction between 0 and 1 '
                    '(usearch only) [default: %default]')),

    make_option('--targetalnfract', type='float', default=0.0,
                help=(
                    'Min percent of the target/reference seq that must match '
                    'to consider a match, expressed as a fraction between 0 '
                    'and 1 (usearch only) [default: %default]')),

    make_option('--max_accepts', type='int', default=1,
                help="max_accepts value (usearch only) [default: %default]"),

    make_option('--max_rejects', type='int', default=32,
                help="max_rejects value to (usearch only) [default: %default]"),

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
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    params = eval(str(opts))
    params['save_uc_files'] = True

    if opts.assignment_method == 'usearch':
        parallel_runner = ParallelDatabaseMapperUsearch(
            cluster_jobs_fp=opts.cluster_jobs_fp,
            jobs_to_start=opts.jobs_to_start,
            retain_temp_files=opts.retain_temp_files,
            suppress_polling=opts.suppress_polling,
            seconds_to_sleep=opts.seconds_to_sleep)
        parallel_runner(opts.input_seqs_filepath,
                        opts.output_dir,
                        params,
                        job_prefix=opts.job_prefix,
                        poll_directly=opts.poll_directly,
                        suppress_submit_jobs=opts.suppress_submit_jobs)
    elif opts.assignment_method == 'blat':
        parallel_runner = ParallelDatabaseMapperBlat(
            cluster_jobs_fp=opts.cluster_jobs_fp,
            jobs_to_start=opts.jobs_to_start,
            retain_temp_files=opts.retain_temp_files,
            suppress_polling=opts.suppress_polling,
            seconds_to_sleep=opts.seconds_to_sleep)
        parallel_runner(opts.input_seqs_filepath,
                        opts.output_dir,
                        params,
                        job_prefix=opts.job_prefix,
                        poll_directly=opts.poll_directly,
                        suppress_submit_jobs=opts.suppress_submit_jobs)
    elif opts.assignment_method == 'bwa-short':
        # cast max_diff to an int if it's being passed as an int
        if params['max_diff'] is not None and params['max_diff'] > 1.0:
            params['max_diff'] = int(params['max_diff'])
        parallel_runner = ParallelDatabaseMapperBwaShort(
            cluster_jobs_fp=opts.cluster_jobs_fp,
            jobs_to_start=opts.jobs_to_start,
            retain_temp_files=opts.retain_temp_files,
            suppress_polling=opts.suppress_polling,
            seconds_to_sleep=opts.seconds_to_sleep)
        parallel_runner(opts.input_seqs_filepath,
                        opts.output_dir,
                        params,
                        job_prefix=opts.job_prefix,
                        poll_directly=opts.poll_directly,
                        suppress_submit_jobs=opts.suppress_submit_jobs)
    else:
        # other -- shouldn't be able to get here as a KeyError would have
        # been raised earlier
        raise ValueError(
            "Unknown read mapping method: %s" %
            opts.assignment_method)


if __name__ == "__main__":
    main()
