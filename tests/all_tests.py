#!/usr/bin/env python
"""Run all tests.
"""
from os import walk, environ
from subprocess import Popen, PIPE, STDOUT
from os.path import join, abspath, dirname, split
from glob import glob
from qiime.util import get_qiime_scripts_dir
import re

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Rob Knight","Greg Caporaso"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

test_dir = abspath(dirname(__file__))

unittest_good_pattern = re.compile('OK\s*$')
application_not_found_pattern = re.compile('ApplicationNotFoundError')
python_name = 'python'
bad_tests = []
missing_application_tests = []

# Run through all of QIIME's unit tests, and keep track of any files which
# fail unit tests.

unittest_names = []

for root, dirs, files in walk(test_dir):
    for name in files:
        if name.startswith('test_') and name.endswith('.py'):
            unittest_names.append(join(root,name))

unittest_names.sort()

for unittest_name in unittest_names:
    print "Testing %s:\n" % unittest_name
    command = '%s %s -v' % (python_name, unittest_name)
    result = Popen(command,shell=True,universal_newlines=True,\
                   stdout=PIPE,stderr=STDOUT).stdout.read()
    print result
    if not unittest_good_pattern.search(result):
        if application_not_found_pattern.search(result):
            missing_application_tests.append(unittest_name)
        else:
            bad_tests.append(unittest_name)


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
    bad_scripts = []

    for script_name in script_names:
        script_good_pattern = re.compile('^Usage: %s' % split(script_name)[1])
        print "Testing %s." % script_name
        command = '%s %s -h' % (python_name, script_name)
        result = Popen(command,shell=True,universal_newlines=True,\
                       stdout=PIPE,stderr=STDOUT).stdout.read()
        if not script_good_pattern.search(result):
            bad_scripts.append(script_name)

if bad_tests:
    print "\nFailed the following unit tests.\n%s" % '\n'.join(bad_tests)
    
if missing_application_tests:
    print "\nFailed the following unit tests, in part or whole due "+\
    "to missing external applications.\nDepending on the QIIME features "+\
    "you plan to use, this may not be critical.\n%s"\
     % '\n'.join(missing_application_tests)
     
if not script_directory_found:
        print "\nCritical error: Failed to test scripts because the script directory could not be found.\n The most likely explanation for this failure is that you've installed QIIME using setup.py, and forgot to specify the qiime_scripts_dir in your qiime_config file. This value shoud be set either to the directory you provided for --install-scripts, or /usr/local/bin if no value was provided to --install-scripts."
else:
    if bad_scripts:
        print "\nFailed the following script tests.\n%s" % '\n'.join(bad_scripts)
     
    if not (bad_tests or missing_application_tests or bad_scripts):
        print "\nAll tests passed successfully."
