#!/usr/bin/env python
"""Run all tests.
"""
from os import walk, popen3, environ
from os.path import join
import re

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

good_pattern = re.compile('OK\s*$') 
start_dir = '.'
python_name = 'python'

bad_tests = []
filenames = []

for root, dirs, files in walk(start_dir):
    for name in files:
        if name.startswith('test_') and name.endswith('.py'):
            filenames.append(join(root,name))
filenames.sort()

for filename in filenames:
    print "Testing %s:\n" % filename
    result = popen3('%s %s -v' % (python_name, filename))[2].read()
    print result
    if not good_pattern.search(result):
        bad_tests.append(filename)

if bad_tests:
    print "Failed the following tests:\n%s" % '\n'.join(bad_tests)
else:
    print "All tests passed successfully."
