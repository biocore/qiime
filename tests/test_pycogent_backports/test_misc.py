#!/usr/bin/env python
# File created on 15 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"
 
from os import rmdir
from os.path import exists
from cogent.util.unit_test import TestCase, main
from qiime.pycogent_backports.misc import get_random_directory_name

class UtilsTests(TestCase):
    """Tests of individual functions in utils"""
    
    def test_get_random_directory_name(self):
        """get_random_directory_name functions as expected """
        # repeated calls yield different directory names
        dirs = []
        for i in range(100):
            d = get_random_directory_name(suppress_mkdir=True)
            self.assertTrue(d not in dirs)
            dirs.append(d)
            
        actual = get_random_directory_name(suppress_mkdir=True)
        self.assertFalse(exists(actual),'Random dir exists: %s' % actual)
        self.assertTrue(actual.startswith('/'),\
         'Random dir is not a full path: %s' % actual)
        
        # prefix, suffix and output_dir are used as expected
        actual = get_random_directory_name(suppress_mkdir=True,prefix='blah',\
            output_dir='/tmp/',suffix='stuff')
        self.assertTrue(actual.startswith('/tmp/blah2'),\
         'Random dir does not begin with output_dir + prefix '+\
         '+ 2 (where 2 indicates the millenium in the timestamp): %s' % actual)
        self.assertTrue(actual.endswith('stuff'),\
         'Random dir does not end with suffix: %s' % actual)
        
        # changing rand_length functions as expected
        actual1 = get_random_directory_name(suppress_mkdir=True)
        actual2 = get_random_directory_name(suppress_mkdir=True,\
            rand_length=10)
        actual3 = get_random_directory_name(suppress_mkdir=True,\
            rand_length=0)
        self.assertTrue(len(actual1) > len(actual2) > len(actual3),\
         "rand_length does not affect directory name lengths "+\
         "as expected:\n%s\n%s\n%s" % (actual1,actual2,actual3))
         
        # changing the timestamp pattern functions as expected
        actual1 = get_random_directory_name(suppress_mkdir=True)
        actual2 = get_random_directory_name(suppress_mkdir=True,\
            timestamp_pattern='%Y')
        self.assertNotEqual(actual1,actual2)
        self.assertTrue(len(actual1)>len(actual2),\
            'Changing timestamp_pattern does not affect directory name')
        # empty string as timestamp works
        actual3 = get_random_directory_name(suppress_mkdir=True,\
            timestamp_pattern='')
        self.assertTrue(len(actual2) > len(actual3))
         
        # creating the directory works as expected 
        actual = get_random_directory_name(output_dir='/tmp/',\
         prefix='get_random_directory_test')
        self.assertTrue(exists(actual))
        rmdir(actual)

if __name__ == "__main__":
    main()