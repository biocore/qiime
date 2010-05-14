#!/usr/bin/env python
from cogent.util.unit_test import TestCase, main
from qiime.trim_sff_primers import (get_technical_lengths) 
"""Tests of the trim_sff_primers.py file.

Note: this is presently just a stub that tests import and the one function that
isn't in the main block: this needs to be refactored.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

class TopLevelTests(TestCase):
    """Top-level tests of functions in trim_sff_primers"""

    def setUp(self):
        """Define a few standard variables"""
        self.with_linker  = """#SampleID\tKEY_SEQ\tBARCODE\tLINKER\tPRIMER\n"""
        self.without_linker  = """#SampleID\tKEY_SEQ\tBARCODE\tLx\tPRIMER\n"""
        self.samples = """a\tATGC\tCCC\tC\tCCCC\nb\tATGC\tGG\tAAA\tCC\n"""
        

    def test_get_technical_lengths(self):
        """get_technical_lengths should return correct dict of sample:length"""
        self.assertEqual(get_technical_lengths(
            (self.with_linker+self.samples).splitlines()),
            {'a':12, 'b':11})
        self.assertEqual(get_technical_lengths(
            (self.without_linker+self.samples).splitlines()), {'a':11, 'b':8})

if __name__ == '__main__':
    main()
