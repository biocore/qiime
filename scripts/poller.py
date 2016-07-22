#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import make_option
from qiime.parallel.poller import poller, get_function_handle
from qiime.util import parse_command_line_parameters

script_info = {}
script_info['brief_description'] = """Poller for parallel QIIME scripts."""
script_info[
    'script_description'] = """Script for polling parallel runs to check completion."""
script_info['script_usage'] = [(
    "Poller example",
    "Runs the poller, which checks for the existence of two input files "
    "(file1.txt and file2.txt) and merges their contents. A cleanup file is "
    "provided that instructs the poller to remove the newly merged file.",
    "%prog -f run_complete.txt -m poller_test_completed.txt -d clean_up.txt"
)]
script_info['version'] = __version__
script_info['output_description'] = "No output created."

script_info['required_options'] = [
    make_option('-f', '--check_run_complete_file',
                help='path to file containing a list of files that must exist to '
                'declare a run complete [REQUIRED]')
]

script_info['optional_options'] = [
    make_option('-r', '--check_run_complete_f',
                help='function which returns True when run is completed '
                '[default: %default]',
                default='qiime.parallel.poller.basic_check_run_complete_f'),
    make_option('-p', '--process_run_results_f',
                help='function to be called when runs complete [default: %default]',
                default='qiime.parallel.poller.basic_process_run_results_f'),
    make_option('-m', '--process_run_results_file',
                help='path to file containing a map of tmp filepaths which should'
                ' be written to final output filepaths [default: %default]'),
    make_option('-c', '--clean_up_f',
                help='function called after processing result [default: %default]',
                default='qiime.parallel.poller.basic_clean_up_f'),
    make_option('-d', '--clean_up_file',
                help='List of files and directories to remove after run'
                ' [default: %default]'),
    make_option('-t', '--time_to_sleep', type='int',
                help='time to wait between calls to status_callback_f'
                ' (in seconds) [default: %default]', default=3)
]


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    poller(get_function_handle(opts.check_run_complete_f),
           get_function_handle(opts.process_run_results_f),
           get_function_handle(opts.clean_up_f),
           list(open(opts.check_run_complete_file)),
           list(open(opts.process_run_results_file)),
           list(open(opts.clean_up_file)),
           seconds_to_sleep=opts.time_to_sleep)


if __name__ == "__main__":
    main()
