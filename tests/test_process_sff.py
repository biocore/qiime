#!/usr/bin/env python
import os
import shutil
import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.app.util import ApplicationNotFoundError
from cogent.util.misc import app_path
from qiime.process_sff import (make_flow_txt,make_fna, make_qual, prep_sffs_in_dir,convert_Ti_to_FLX)
from qiime.util import get_qiime_project_dir

"""Tests of the process_sff.py file.
"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
#remember to add yourself if you make changes
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

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
        self.outpath = os.path.join(self.sff_dir,'test')
        shutil.copy(sff_original_fp, self.sff_fp)

    def tearDown(self):
        shutil.rmtree(self.sff_dir)

    def test_convert_Ti_to_FLX(self):
        """test_make_flow_txt should make flowgram file as expected"""
        convert_Ti_to_FLX(self.sff_fp,self.outpath+'_FLX.sff')
        expected_fp = os.path.join(self.sff_dir,'test_FLX.sff')
        fsize=os.path.getsize(expected_fp)
        self.failIfEqual(0, fsize)
    
    def test_make_flow_txt(self):
        """test_make_flow_txt should make flowgram file as expected"""
        make_flow_txt(self.sff_fp,self.outpath)
        expected_fp = os.path.join(self.sff_dir,'test.txt')
        observed = open(expected_fp).read()
        self.assertEqual(observed, flow_expected)

    def test_make_fna(self):
        """test_make_fna should make fasta file as expected"""
        make_fna(self.sff_fp,self.outpath)
        expected_fp = os.path.join(self.sff_dir, 'test.fna')
        observed = open(expected_fp).read()
        expected = (
            '>FA6P1OK01CGMHQ length=48 xy=0892_1356 region=1 run=R_2008_05_28_17_11_38_\n'
            'ATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGC\n'
            )
        self.assertEqual(observed, expected)

    def test_make_qual(self):
        """test_make_qual should make qual file as expected"""
        make_qual(self.sff_fp,self.outpath)
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
        prep_sffs_in_dir(self.sff_dir,False,self.sff_dir,False)
        
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

flow_expected='''Common Header:\n  Magic Number:  0x2E736666\n  Version:       0001\n  Index Offset:  1504\n  Index Length:  706\n  # of Reads:    1\n  Header Length: 440\n  Key Length:    4\n  # of Flows:    400\n  Flowgram Code: 1\n  Flow Chars:    TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\n  Key Sequence:  TCAG\n\n>FA6P1OK01CGMHQ\n  Run Prefix:   R_2008_05_28_17_11_38_\n  Region #:     1\n  XY Location:  0892_1356\n\n  Run Name:       R_2008_05_28_17_11_38_FLX02070135_adminrig_KnightLauber\n  Analysis Name:  /data/2008_05_28/R_2008_05_28_17_11_38_FLX02070135_adminrig_KnightLauber/D_2008_05_28_21_13_06_FLX02070135_KnightLauber_FullAnalysisAmplicons\n  Full Path:      /data/2008_05_28/R_2008_05_28_17_11_38_FLX02070135_adminrig_KnightLauber/D_2008_05_28_21_13_06_FLX02070135_KnightLauber_FullAnalysisAmplicons/../D_2008_05_29_13_52_01_FLX02070135_Knight_Lauber_jones_SignalProcessingAmplicons\n\n  Read Header Len:  32\n  Name Length:      14\n  # of Bases:       77\n  Clip Qual Left:   5\n  Clip Qual Right:  52\n  Clip Adap Left:   0\n  Clip Adap Right:  0\n\nFlowgram:\t1.04\t0.00\t1.01\t0.00\t0.00\t0.96\t0.00\t1.02\t0.00\t1.02\t0.00\t0.00\t0.99\t0.00\t1.00\t0.00\t1.00\t0.00\t0.00\t1.00\t0.00\t1.10\t0.00\t1.08\t0.00\t0.00\t1.46\t0.00\t0.88\t0.18\t0.00\t2.69\t1.01\t0.08\t0.96\t0.00\t0.02\t0.92\t0.08\t0.00\t0.98\t0.68\t0.00\t0.89\t0.00\t0.00\t1.15\t0.00\t1.13\t0.00\t0.02\t1.12\t0.05\t0.15\t1.84\t0.00\t1.10\t0.00\t2.47\t0.96\t0.86\t1.06\t0.00\t1.96\t0.12\t0.93\t0.13\t1.65\t1.06\t0.06\t0.00\t0.99\t0.00\t0.00\t1.87\t0.44\t1.08\t0.00\t3.25\t0.09\t0.97\t0.50\t1.00\t1.72\t0.07\t0.00\t0.92\t0.58\t0.00\t0.00\t0.59\t0.06\t0.11\t0.09\t0.07\t0.06\t0.16\t0.00\t0.24\t0.03\t0.00\t0.12\t0.06\t0.16\t0.00\t0.18\t0.00\t0.00\t0.14\t0.00\t0.15\t0.00\t0.18\t0.00\t0.03\t0.14\t0.03\t0.13\t0.01\t0.19\t0.00\t0.02\t0.33\t0.05\t0.00\t0.16\t0.10\t0.35\t0.01\t0.21\t0.04\t0.09\t0.18\t0.13\t0.19\t0.00\t0.10\t0.51\t0.26\t0.00\t0.23\t0.19\t0.27\t0.01\t0.29\t0.05\t0.14\t0.17\t0.16\t0.18\t0.27\t0.09\t0.26\t0.10\t0.18\t0.23\t0.15\t0.22\t0.13\t0.37\t0.11\t0.11\t0.26\t0.59\t0.14\t0.06\t0.33\t0.34\t0.26\t0.05\t0.27\t0.44\t0.19\t0.10\t0.35\t0.27\t0.15\t0.34\t0.28\t0.45\t0.14\t0.16\t0.34\t0.27\t0.12\t0.07\t0.25\t0.18\t0.12\t0.04\t0.23\t0.16\t0.12\t0.05\t0.20\t0.16\t0.11\t0.03\t0.21\t0.16\t0.10\t0.02\t0.21\t0.16\t0.12\t0.02\t0.20\t0.15\t0.10\t0.02\t0.23\t0.15\t0.11\t0.02\t0.22\t0.14\t0.09\t0.02\t0.20\t0.13\t0.09\t0.01\t0.19\t0.13\t0.08\t0.02\t0.17\t0.12\t0.08\t0.03\t0.17\t0.09\t0.08\t0.01\t0.14\t0.09\t0.07\t0.01\t0.15\t0.09\t0.06\t0.01\t0.13\t0.08\t0.06\t0.00\t0.13\t0.08\t0.05\t0.02\t0.12\t0.07\t0.05\t0.01\t0.11\t0.07\t0.05\t0.00\t0.10\t0.07\t0.05\t0.01\t0.11\t0.08\t0.04\t0.00\t0.10\t0.06\t0.05\t0.01\t0.09\t0.06\t0.04\t0.01\t0.08\t0.07\t0.05\t0.00\t0.08\t0.06\t0.05\t0.00\t0.09\t0.06\t0.04\t0.00\t0.09\t0.06\t0.04\t0.01\t0.08\t0.06\t0.04\t0.00\t0.09\t0.06\t0.03\t0.00\t0.09\t0.06\t0.02\t0.00\t0.09\t0.06\t0.04\t0.00\t0.08\t0.05\t0.03\t0.00\t0.07\t0.05\t0.02\t0.00\t0.08\t0.04\t0.03\t0.00\t0.07\t0.04\t0.03\t0.00\t0.07\t0.05\t0.02\t0.00\t0.07\t0.05\t0.02\t0.00\t0.06\t0.04\t0.02\t0.00\t0.06\t0.03\t0.03\t0.00\t0.08\t0.02\t0.00\t0.00\t0.07\t0.03\t0.01\t0.00\t0.06\t0.03\t0.02\t0.00\t0.05\t0.03\t0.02\t0.00\t0.05\t0.03\t0.01\t0.00\t0.06\t0.02\t0.00\t0.00\t0.05\t0.01\t0.01\t0.00\t0.04\t0.01\t0.01\t0.00\t0.04\t0.01\t0.01\t0.00\t0.05\t0.01\t0.00\t0.00\t0.04\t0.02\t0.01\t0.00\t0.03\t0.02\t0.01\t0.00\t0.03\t0.01\t0.00\t0.00\t0.03\t0.00\t0.00\t0.00\t0.03\t0.00\t0.00\t0.00\t0.02\t0.00\nFlow Indexes:\t1\t3\t6\t8\t10\t13\t15\t17\t20\t22\t24\t27\t29\t32\t32\t32\t33\t35\t38\t41\t42\t44\t47\t49\t52\t55\t55\t57\t59\t59\t60\t61\t62\t64\t64\t66\t68\t68\t69\t72\t75\t75\t77\t79\t79\t79\t81\t82\t83\t84\t84\t87\t88\t91\t102\t126\t130\t138\t140\t145\t153\t157\t161\t164\t166\t171\t175\t179\t183\t187\t191\t195\t199\t203\t211\t215\t219\nBases:\ttcagATCTGAGCTGGGTCATAGCTGCCTCCGTAGGAGGTGCCTCCCTACGGCgcnnnannnnngnnnnnnnnnnnnn\nQuality Scores:\t32\t32\t32\t32\t32\t32\t32\t32\t32\t32\t32\t25\t25\t21\t21\t21\t28\t32\t32\t31\t30\t30\t32\t32\t32\t33\t31\t25\t18\t18\t20\t18\t32\t30\t28\t23\t22\t22\t24\t28\t18\t19\t18\t16\t16\t16\t17\t18\t13\t17\t27\t21\t20\t21\t0\t0\t0\t17\t0\t0\t0\t0\t0\t17\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n'''


if __name__ == '__main__':
    main()

