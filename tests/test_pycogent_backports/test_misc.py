#!/usr/bin/env python

from os import rmdir
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files, get_random_directory_name
from qiime.pycogent_backports.misc import (parse_command_line_parameters, 
                                           handle_error_codes, 
                                           create_dir)
from numpy import array

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2009, The Cogent Project"
__credits__ = ["Rob Knight", "Amanda Birmingham", "Sandra Smit",
                    "Zongzhi Liu", "Peter Maxwell", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0.dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class CommandLineParserTests(TestCase):
    
    
    def setUp(self):
        """ """
        self.files_to_remove = []
        
    def tearDown(self):
        """ """
        map(rmdir,self.files_to_remove)
    
    
    def test_create_dir(self):
        """create_dir creates dir and fails meaningful."""

        tmp_dir_path = get_random_directory_name()
        tmp_dir_path2 = get_random_directory_name(suppress_mkdir=True)
        tmp_dir_path3 = get_random_directory_name(suppress_mkdir=True)

        self.files_to_remove.append(tmp_dir_path)
        self.files_to_remove.append(tmp_dir_path2)
        self.files_to_remove.append(tmp_dir_path3)

        # create on existing dir raises OSError if fail_on_exist=True
        self.assertRaises(OSError, create_dir, tmp_dir_path,
                          fail_on_exist=True)
        self.assertEquals(create_dir(tmp_dir_path,
                                     fail_on_exist=True,
                                     handle_errors_externally=True), 1)

        # return should be 1 if dir exist and fail_on_exist=False 
        self.assertEqual(create_dir(tmp_dir_path, fail_on_exist=False), 1)

        # if dir not there make it and return always 0
        self.assertEqual(create_dir(tmp_dir_path2), 0)
        self.assertEqual(create_dir(tmp_dir_path3, fail_on_exist=True), 0)

    def test_handle_error_codes(self):
        """handle_error_codes raises the right error."""

        self.assertRaises(OSError, handle_error_codes, "test", False,1)
        self.assertEqual(handle_error_codes("test", True, 1), 1)
        self.assertEqual(handle_error_codes("test", False, 0), 0)
        self.assertEqual(handle_error_codes("test"), 0)
    
    def test_parse_command_line_parameters(self):
        """parse_command_line_parameters returns without error
        
            There is not a lot of detailed testing that can be done here,
            so the basic functionality is tested.
        """
        option_parser, opts, args = parse_command_line_parameters(
          script_description="My script",
          script_usage=[('Print help','%prog -h','')],
          version='1.0',help_on_no_arguments=False)
        self.assertEqual(len(args),0)
        
        d = {'script_description':"My script",\
             'script_usage':[('Print help','%prog -h','')],\
             'version':'1.0',
             'help_on_no_arguments':False}
        option_parser, opts, args = parse_command_line_parameters(**d)
        self.assertEqual(len(args),0)

#run tests on command-line invocation

if __name__ == '__main__':
    main()
