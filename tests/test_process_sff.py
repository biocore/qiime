#!/usr/bin/env python
import os
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationNotFoundError
from cogent.util.misc import app_path
from qiime.process_sff import (make_fna, make_qual, prep_sffs_in_dir)
from qiime.util import get_qiime_project_dir

"""Tests of the process_sff.py file.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"

class TopLevelTests(TestCase):
    """Top-level tests of functions in process_sff"""
    
    def setUp(self):
        """Create temporary directory of SFF files"""

        # Cannot use get_qiime_project_dir() due to test errors in virtual box
        test_dir = os.path.dirname(os.path.abspath(__file__))
        sff_original_fp =  os.path.join(test_dir, 'sra_test_files', 'test.sff')

        # copy sff file to working directory
        self.sff_dir = tempfile.mkdtemp()
        self.sff_fp = os.path.join(self.sff_dir, 'test.sff')
        shutil.copy(sff_original_fp, self.sff_fp)

    def tearDown(self):
        shutil.rmtree(self.sff_dir)

    def test_make_fna(self):
        """test_make_fna should make fasta file as expected"""
        make_fna(self.sff_fp)
        expected_fp = os.path.join(self.sff_dir, 'test.fna')
        observed = open(expected_fp).read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n'
            'ATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n'
            )
        self.assertEqual(observed, expected)

    def test_make_qual(self):
        """test_make_qual should make qual file as expected"""
        make_qual(self.sff_fp)
        expected_fp = os.path.join(self.sff_dir, 'test.qual')
        observed = open(expected_fp).read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n'
            '32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 '
            '20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21\n'
            )
        self.assertEqual(observed, expected)

    def test_prep_sffs_in_dir(self):
        """test_prep_sffs_in_dir should make fasta/qual from sffs."""
        prep_sffs_in_dir(self.sff_dir)

        # Check fna file
        expected_fp = os.path.join(self.sff_dir, 'test.fna')
        observed = open(expected_fp).read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n'
            'ATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n'
            )
        self.assertEqual(observed, expected)

        # Check qual file
        expected_fp = os.path.join(self.sff_dir, 'test.qual')
        observed = open(expected_fp).read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n'
            '32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 '
            '20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21\n'
            )
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    main()
