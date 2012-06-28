#!/usr/bin/env python
"""Run all tests.
"""
from os import walk, environ
from subprocess import Popen, PIPE, STDOUT
from os.path import join, abspath, dirname, split
from glob import glob
import re
from sys import exit
from cogent.app.util import get_tmp_filename
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                       load_qiime_config,qiime_system_call,get_qiime_scripts_dir,
                       make_option)
from qiime.test import run_script_usage_tests

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Rob Knight","Greg Caporaso", "Jai Ram Rideout"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = []
script_info['optional_options'] = [
 make_option('--suppress_unit_tests',
             action='store_true',
             help='suppress unit tests [default: %default]',
             default=False),
 make_option('--suppress_script_tests',
             action='store_true',
             help='suppress script tests [default: %default]',
             default=False),
 make_option('--suppress_script_usage_tests',
             action='store_true',
             help='suppress script usage tests [default: %default]',
             default=False),
]
script_info['version'] = __version__
script_info['help_on_no_arguments'] = False

qiime_config = load_qiime_config()

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    if (opts.suppress_unit_tests and \
       opts.suppress_script_tests and \
       opts.suppress_script_usage_tests):
       option_parser.error("You're suppressing all three test types. Nothing to run.")

    test_dir = abspath(dirname(__file__))

    unittest_good_pattern = re.compile('OK\s*$')
    application_not_found_pattern = re.compile('ApplicationNotFoundError')
    python_name = 'python'
    bad_tests = []
    missing_application_tests = []

    # Run through all of QIIME's unit tests, and keep track of any files which
    # fail unit tests.
    if not opts.suppress_unit_tests:
        unittest_names = []

        for root, dirs, files in walk(test_dir):
            for name in files:
                if name.startswith('test_') and name.endswith('.py'):
                    unittest_names.append(join(root,name))

        unittest_names.sort()

        for unittest_name in unittest_names:
            print "Testing %s:\n" % unittest_name
            command = '%s %s -v' % (python_name, unittest_name)
            stdout, stderr, return_value = qiime_system_call(command)
            print stderr
            if not unittest_good_pattern.search(stderr):
                if application_not_found_pattern.search(stderr):
                    missing_application_tests.append(unittest_name)
                else:
                    bad_tests.append(unittest_name)


    bad_scripts = []
    if not opts.suppress_script_tests:
        # Run through all of QIIME's scripts, and pass -h to each one. If the
        # resulting stdout does not being with the Usage text, that is an 
        # indicator of something being wrong with the script. Issues that would
        # cause that are bad import statements in the script, SyntaxErrors, or 
        # other failures prior to running qiime.util.parse_command_line_parameters.

        try:
            scripts_dir = get_qiime_scripts_dir()
            script_directory_found = True
        except AssertionError:
            script_directory_found = False


        if script_directory_found:
            script_names = []
            script_names = glob('%s/*py' % scripts_dir)
            script_names.sort()

            for script_name in script_names:
                script_good_pattern = re.compile('^Usage: %s' % split(script_name)[1])
                print "Testing %s." % script_name
                command = '%s %s -h' % (python_name, script_name)
                stdout, stderr, return_value = qiime_system_call(command)
                if not script_good_pattern.search(stdout):
                    bad_scripts.append(script_name)
    
    num_script_usage_example_failures = 0
    qiime_test_data_dir = qiime_config['qiime_test_data_dir']
    if not opts.suppress_script_usage_tests and qiime_test_data_dir != None:
        # Run the script usage testing functionality
        script_usage_result_summary, num_script_usage_example_failures = \
         run_script_usage_tests(
               qiime_test_data_dir=qiime_test_data_dir,
               qiime_scripts_dir=qiime_config['qiime_scripts_dir'],
               working_dir=qiime_config['temp_dir'],
               verbose=True,
               tests=None, # runs all
               failure_log_fp=None,
               force_overwrite=True)

    print "==============\nResult summary\n=============="

    if not opts.suppress_unit_tests:
        print "\nUnit test result summary\n------------------------\n"
        if bad_tests:
            print "\nFailed the following unit tests.\n%s" % '\n'.join(bad_tests)
    
        if missing_application_tests:
            print "\nFailed the following unit tests, in part or whole due "+\
            "to missing external applications.\nDepending on the QIIME features "+\
            "you plan to use, this may not be critical.\n%s"\
             % '\n'.join(missing_application_tests)
        
        if not (missing_application_tests or bad_tests):
            print "\nAll unit tests passed.\n\n"
     
    if not opts.suppress_script_tests:
        print "\nBasic script test result summary\n--------------------------------\n"
        if not script_directory_found:
            print "Critical error: Failed to test scripts because the script directory could not be found.\n The most likely explanation for this failure is that you've installed QIIME using setup.py, and forgot to specify the qiime_scripts_dir in your qiime_config file. This value shoud be set either to the directory you provided for --install-scripts, or /usr/local/bin if no value was provided to --install-scripts."
        else:
            if bad_scripts:
                print "Failed the following basic script tests.\n%s" % '\n'.join(bad_scripts)
            else:
                print "All basic script tests passed successfully.\n"
    
    qiime_test_data_dir_exists = True
    if not opts.suppress_script_usage_tests:
        if qiime_test_data_dir:
            print "\nScript usage test result summary\n------------------------------------\n"
            print script_usage_result_summary
        else:
            print "\nCould not run script usage tests because qiime_test_data_dir is not defined in your qiime_config."
            qiime_test_data_dir_exists = False
        print ""

    # If any of the unit tests, script tests, or script usage tests fail, or if
    # we have any missing application errors or a missing QIIME test data dir
    # if script usage tests weren't suppressed, use return code 1 (as python's
    # unittest module does to indicate one or more failures).
    return_code = 1
    if (len(bad_tests) == 0 and len(missing_application_tests) == 0 and
        len(bad_scripts) == 0 and num_script_usage_example_failures == 0 and
        qiime_test_data_dir_exists):
        return_code = 0
    return return_code


if __name__ == "__main__":
    exit(main())
