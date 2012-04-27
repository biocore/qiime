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
 

from os.path import isdir, split, join, abspath
from os import chdir, getcwd
from shutil import copytree, rmtree
from glob import glob
from site import addsitedir
from sys import path
from cogent.util.misc import remove_files
from qiime.util import (parse_command_line_parameters,
                        make_option,
                        load_qiime_config,
                        qiime_system_call)


qiime_config = load_qiime_config()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
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
 make_option('-w','--working_dir',default=qiime_config['temp_dir'],
             help='the tests to run [default: %default]',type='existing_dirpath'),
 make_option('-q','--qiime_scripts_dir',default=qiime_config['qiime_scripts_dir'],
             help='directory containing scripts to test [default: %default]',
             type='existing_dirpath'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    tests = opts.tests
    qiime_test_data_dir = abspath(opts.qiime_test_data_dir)
    qiime_scripts_dir = opts.qiime_scripts_dir
    working_dir = join(opts.working_dir,'script_usage_tests')
    verbose = opts.verbose
    failure_log_fp = abspath(opts.failure_log_fp)
    
    if tests == None:
        tests = [split(d)[1] for d in glob('%s/*' % qiime_test_data_dir) if isdir(d)]
    else:
        tests = tests.split(',')
    
    if verbose:
        print 'Tests to run:\n %s' % ' '.join(tests)
    
    addsitedir(qiime_scripts_dir)
    
    failed_tests = []
    warnings = []
    total_tests = 0
    for test in tests:
        
        # import the usage examples - this is possible because we added 
        # qiime_scripts_dir to the PYTHONPATH above
        script_fn = '%s/%s.py' % (qiime_scripts_dir,test)
        script = __import__(test)
        usage_examples = script.script_info['script_usage']
        
        if verbose:
            print 'Testing %d usage examples from: %s.py' % (len(usage_examples),script_fn)
        
        # init the test environment
        test_input_dir = '%s/%s' % (qiime_test_data_dir,test)
        test_working_dir = '%s/%s' % (working_dir,test)
        copytree(test_input_dir,test_working_dir)
        chdir(test_working_dir)
        
        # remove pre-exisitng output files if any
        try:
            script_usage_output_to_remove = script.script_info['script_usage_output_to_remove']
        except KeyError:
            script_usage_output_to_remove = []
        for e in script_usage_output_to_remove:
            rmtree(e.replace('$PWD',getcwd()),ignore_errors=True)
        
        if verbose:
            print ' Running tests in: %s' % getcwd()
            print ' Tests:'
        
        for usage_example in usage_examples:
            if '%prog' not in usage_example[2]:
                warnings.append('%s usage examples do not all use %%prog to represent the command name. You may not be running the version of the command that you think you are!' % test)
            cmd = usage_example[2].replace('%prog',script_fn)
            if verbose:
                print '  %s' % cmd
            stdout, stderr, return_value = qiime_system_call(cmd)
            total_tests += 1
            if return_value != 0:
                failed_tests.append((cmd, stdout, stderr, return_value))
        
        if verbose:
            print ''
    
    failure_log_f = open(failure_log_fp,'w')
    if len(failed_tests) == 0:
        failure_log_f.write('All tests passed.')
    else:
        i = 0
        for cmd, stdout, stderr, return_value in failed_tests:
            failure_log_f.write('**Failed test %d:\n%s\n\nReturn value: %d\n\nStdout:\n%s\n\nStderr:\n%s\n\n' % (i,cmd,return_value, stdout, stderr))
    failure_log_f.close()
    
    
    if warnings:
        print 'Warnings:'
        for warning in warnings:
            print ' ' + warning
        print ''
    print 'Ran %d commands to test %d scripts. %d of these commands failed. Failures are summarized in %s.' % (total_tests,len(tests),len(failed_tests),failure_log_fp)
    
    rmtree(working_dir)



if __name__ == "__main__":
    main()