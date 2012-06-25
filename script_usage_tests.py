#!/usr/bin/env python
# File created on 15 Mar 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 


from qiime.util import (parse_command_line_parameters,
                        make_option, load_qiime_config)
from qiime.test import run_script_usage_tests


qiime_config = load_qiime_config()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("Run a subset of the interface tests in verbose mode","Run interface tests for the add_taxa.py and make_otu_table.py scripts. This illustrates how to run from the qiime_test_dir directory.","%prog -i $PWD/ -l $HOME/qime_script_tests.log -t add_taxa,make_otu_table -v"))
script_info['script_usage'].append(("Run all of the interface tests","Run all script interface tests.  This illustrates how to run from the qiime_test_dir directory.","%prog -i $PWD/ -l $HOME/all_qime_script_tests.log"))
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--qiime_test_data_dir',type="existing_dirpath",
             help='the directory containing input for script usage examples'),
 make_option('-l','--failure_log_fp',type="new_filepath",
             help='log file to store record of failures'),
]

script_info['optional_options'] = [\
 make_option('-t','--tests',
             help='comma-separated list of the tests to run [default: all]'),
 make_option('-w','--working_dir',default=qiime_config['temp_dir'] or '/tmp/',
             help='directory where the tests should be run [default: %default]',
             type='existing_dirpath'),
 make_option('-q','--qiime_scripts_dir',default=qiime_config['qiime_scripts_dir'],
             help='directory containing scripts to test [default: %default]',
             type='existing_dirpath'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    qiime_test_data_dir = opts.qiime_test_data_dir
    qiime_scripts_dir = opts.qiime_scripts_dir
    working_dir = opts.working_dir
    verbose = opts.verbose
    tests = opts.tests
    if tests != None:
        tests = tests.split(',')
    failure_log_fp = opts.failure_log_fp

    result_summary = run_script_usage_tests(
                           qiime_test_data_dir,
                           qiime_scripts_dir,
                           working_dir,
                           verbose=verbose,
                           tests=tests,
                           failure_log_fp=failure_log_fp)
    if verbose:
        print result_summary



if __name__ == "__main__":
    main()