#!/usr/bin/env python
from cogent.util.unit_test import TestCase, main
from qiime.sra_spreadsheet_to_map_files import (strip_quotes, 
    collect_study_groups, remap_lines, write_map_files) 
"""Tests of the sra_spreadsheet_to_map_files.py file.
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
    """Top-level tests of functions in sra_spreadsheet_to_map_files.py"""

    def test_strip_quotes(self):
        """strip_quotes should strip terminal quotes from field."""
        self.assertEqual(strip_quotes('abc'), 'abc')
        self.assertEqual(strip_quotes('"abc"'), 'abc')

    def test_remap_lines(self):
        """remap_lines should fix some issues with input lines."""
        lines = """POOL_MEMBER_NAME\tBARCODE\tPRIMER\tLINKER\tabc
x\tAA\tGGG\tCC\tx_x
y\tAC\tCCC\tAA\ty_y""".splitlines()
        result = remap_lines(lines[0].split('\t'), 
            [i.split('\t') for i in lines[1:]])
        self.assertEqual(result,
        [['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence'] + 
            'POOL_MEMBER_NAME\tBARCODE\tPRIMER\tLINKER\tabc'.split('\t') +
            ['Description'],
            ['x','AA','CCGGG','x','AA','GGG','CC','x_x','None'],
            ['y','AC','AACCC','y','AC','CCC','AA','y_y','None']])
        

if __name__ == '__main__':
    main()
