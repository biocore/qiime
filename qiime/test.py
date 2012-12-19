#!/usr/bin/env python
# File created on 27 Apr 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

import signal
from os.path import isdir, split, join, abspath, exists
from os import chdir, getcwd
from shutil import copytree, rmtree
from glob import glob
from site import addsitedir
from sys import path
from cogent.util.misc import remove_files
from qiime.util import qiime_system_call


### Code for timing out tests that exceed a time limit
## The test case timing code included in this file is adapted from
## recipes provided at:
##  http://code.activestate.com/recipes/534115-function-timeout/
##  http://stackoverflow.com/questions/492519/timeout-on-a-python-function-call

# to use this, call initiate_timeout(allowed_seconds_per_test) in 
# TestCase.setUp() and then disable_timeout() in TestCase.tearDown()

class TimeExceededError(Exception):
    pass

def initiate_timeout(seconds=60):
    
    def timeout(signum, frame):
        raise TimeExceededError,\
         "Test failed to run in allowed time (%d seconds)." % seconds
    
    signal.signal(signal.SIGALRM, timeout)
    # set the 'alarm' to go off in seconds seconds
    signal.alarm(seconds)

def disable_timeout():
    # turn off the alarm
    signal.alarm(0)

### End code for timing out tests that exceed a time limit

class FakeFile(object):
    " A class to convert a string into a file-like object"
    def __init__(self,d=""):
        self.s = d.split('\n')
    
    def __iter__(self):
        return iter(self.s)
    
    def write(self,s):
        self.s += s
    def close(self):
        pass
    def seek(self,seek_position):
        pass


def run_script_usage_tests(qiime_test_data_dir,
                           qiime_scripts_dir,
                           working_dir,
                           verbose=False,
                           tests=None,
                           failure_log_fp=None,
                           force_overwrite=False):
    """ Test script_usage examples when test data is present in qiime_test_data_dir
    
        These tests are currently used with the qiime_test_data repository, which can
        be found at https://github.com/qiime-dev/qiime_test_data

        Returns a result summary string and the number of script usage
        examples (i.e. commands) that failed.
    """
    # process input filepaths and directories
    qiime_test_data_dir = abspath(qiime_test_data_dir)
    working_dir = join(working_dir,'script_usage_tests')
    if force_overwrite and exists(working_dir):
        rmtree(working_dir)
    if failure_log_fp != None:
        failure_log_fp = abspath(failure_log_fp)

    if tests == None:
        tests = [split(d)[1] for d in glob('%s/*' % qiime_test_data_dir) if isdir(d)]
    
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
            print 'Testing %d usage examples from: %s' % (len(usage_examples),script_fn)
        
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
            remove_files([e.replace('$PWD',getcwd())],error_on_missing=False)
        
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
            
    if failure_log_fp:
        failure_log_f = open(failure_log_fp,'w')
        if len(failed_tests) == 0:
            failure_log_f.write('All script interface tests passed.\n')
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
    
    result_summary = 'Ran %d commands to test %d scripts. %d of these commands failed.' % (total_tests,len(tests),len(failed_tests))
    if len(failed_tests) > 0:
        failed_scripts = set([split(e[0].split()[0])[1] for e in failed_tests])
        result_summary += '\nFailed scripts were: %s' % " ".join(failed_scripts)
    if failure_log_fp:
        result_summary += "\nFailures are summarized in %s" % failure_log_fp
    
    rmtree(working_dir)
    
    return result_summary, len(failed_tests)
