#!/usr/bin/env python
from cogent.util.unit_test import TestCase, main
from qiime.process_sff import (make_fna, make_qual, prep_sffs_in_dir) 
from os import remove
"""Tests of the process_sff.py file.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

class TopLevelTests(TestCase):
    """Top-level tests of functions in process_sff"""

    def test_make_fna(self):
        """test_make_fna should make fasta file as expected"""
        make_fna('sra_test_files/test.sff')
        result = open('sra_test_files/test.fna').read()
        self.assertEqual(result, '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\nATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n')
        remove('sra_test_files/test.fna')

    def test_make_qual(self):
        """test_make_qual should make qual file as expected"""
        make_qual('sra_test_files/test.sff')
        result = open('sra_test_files/test.qual').read()
        self.assertEqual(result, '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21\n')
        remove('sra_test_files/test.qual')

    def test_prep_sffs_in_dir(self):
        """test_prep_sffs_in_dir should make fasta/qual from sffs."""
        prep_sffs_in_dir('sra_test_files')
        result = open('sra_test_files/test.qual').read()
        self.assertEqual(result, '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21\n')
        remove('sra_test_files/test.qual')
        result = open('sra_test_files/test.fna').read()
        self.assertEqual(result, '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\nATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n')
        remove('sra_test_files/test.fna')

if __name__ == '__main__':
    main()
