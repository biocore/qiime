#!/usr/bin/env python
# File created on 15 Mar 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 


from qiime.util import (parse_command_line_parameters,
                        make_option, load_qiime_config, get_qiime_temp_dir)
from qiime.test import run_script_usage_tests
from qiime.workflow import generate_log_fp
from os.path import join

qiime_config = load_qiime_config()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("Run a subset of the interface tests in verbose mode","Run interface tests for the count_seqs.py and make_otu_table.py scripts. This illustrates how to run from the qiime_test_dir directory.","%prog -t count_seqs,make_otu_table -v"))
script_info['script_usage'].append(("Run all of the interface tests","Run all script interface tests.  This illustrates how to run from the qiime_test_dir directory.","%prog"))
script_info['output_description']= ""
script_info['required_options'] = []

log_fp_prefix = 'script_test_log'
log_fp_suffix = 'txt'
default_log_fp = generate_log_fp(get_qiime_temp_dir(),
                    basefile_name=log_fp_prefix,
                    suffix=log_fp_suffix,
                    timestamp_pattern='%Y%m%d%H%M%S')
default_log_fp_help_str = join(get_qiime_temp_dir(),
                               '%s_TIMESTAMP.%s' % (log_fp_prefix,log_fp_suffix))

script_info['optional_options'] = [\
 make_option('-t','--tests',
             help='comma-separated list of the tests to run [default: all]'),
 make_option('-w','--working_dir',default=get_qiime_temp_dir(),
             help='directory where the tests should be run [default: %default]',
             type='existing_dirpath'),
 make_option('-l','--failure_log_fp',type="new_filepath",default=default_log_fp,
             help='log file to store record of failures [default: %s]' % default_log_fp_help_str)
]
script_info['version'] = __version__

default_qiime_test_data_dir = qiime_config['qiime_test_data_dir']
if default_qiime_test_data_dir != None:
    script_info['optional_options'].append(
     make_option('-i','--qiime_test_data_dir',type="existing_dirpath",
                 default=default_qiime_test_data_dir,
                 help='the directory containing input for script usage'+\
                 ' examples [default: %default]'))
else:
    script_info['required_options'].append(
     make_option('-i','--qiime_test_data_dir',type="existing_dirpath",
                 help='the directory containing input for script usage examples'))

default_qiime_scripts_dir = qiime_config['qiime_scripts_dir']
if default_qiime_scripts_dir != None:
    script_info['optional_options'].append(
     make_option('-q','--qiime_scripts_dir',
             default=default_qiime_scripts_dir,
             help='directory containing scripts to test [default: %default]',
             type='existing_dirpath'))
else:
    script_info['required_options'].append(
     make_option('-q','--qiime_scripts_dir',
             help='directory containing scripts to test',
             type='existing_dirpath'))

script_info['help_on_no_arguments'] = False

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    qiime_test_data_dir = opts.qiime_test_data_dir
    qiime_scripts_dir = opts.qiime_scripts_dir
    working_dir = opts.working_dir
    verbose = opts.verbose
    tests = opts.tests
    if tests != None:
        tests = [e.rstrip('/') for e in tests.split(',')]
    failure_log_fp = opts.failure_log_fp
    
    result_summary, num_failures = run_script_usage_tests(
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
