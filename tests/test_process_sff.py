#!/usr/bin/env python
import os
from os import remove, getcwd
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
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

class TopLevelTests(TestCase):
    """Top-level tests of functions in process_sff"""
    
    def setUp(self):
        """ """
        qiime_dir = get_qiime_project_dir()
        self.sra_test_files_dir = os.path.join(
            qiime_dir, 'tests', 'sra_test_files')

    def test_make_fna(self):
        """test_make_fna should make fasta file as expected"""
        if not app_path('sffinfo'):
            raise ApplicationNotFoundError, \
             "Can't find sffinfo. Is it installed? Is it in your $PATH?"
        make_fna('%s/test.sff' % self.sra_test_files_dir)
        result = open('%s/test.fna' % self.sra_test_files_dir).read()
        self.assertEqual(result, \
         '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\nATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n')
        remove('%s/test.fna' % self.sra_test_files_dir)

    def test_make_qual(self):
        """test_make_qual should make qual file as expected"""
        if not app_path('sffinfo'):
            raise ApplicationNotFoundError, \
             "Can't find sffinfo. Is it installed? Is it in your $PATH?"
        make_qual('%s/test.sff' % self.sra_test_files_dir)
        result = open('%s/test.qual'% self.sra_test_files_dir).read()
        self.assertEqual(result,\
        '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21\n')
        remove('%s/test.qual' % self.sra_test_files_dir)

    def test_prep_sffs_in_dir(self):
        """test_prep_sffs_in_dir should make fasta/qual from sffs."""
        if not app_path('sffinfo'):
            raise ApplicationNotFoundError, \
             "Can't find sffinfo. Is it installed? Is it in your $PATH?"
        prep_sffs_in_dir(self.sra_test_files_dir)
        result = open('%s/test.qual' % self.sra_test_files_dir).read()
        self.assertEqual(result, \
        '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n32 32 32 32 32 32 32 25 25 21 21 21 28 32 32 31 30 30 32 32 32 33 31 25 18 18 20 18 32 30 28 23 22 22 24 28 18 19 18 16 16 16 17 18 13 17 27 21\n')
        remove('%s/test.qual' % self.sra_test_files_dir)
        result = open('%s/test.fna' % self.sra_test_files_dir).read()
        self.assertEqual(result, '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\nATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n')
        remove('%s/test.fna' % self.sra_test_files_dir)

if __name__ == '__main__':
    main()
