#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from qiime.pycogent_backports.misc import parse_command_line_parameters
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
