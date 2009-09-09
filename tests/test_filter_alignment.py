#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

from cogent.util.unit_test import TestCase, main
from cogent import LoadSeqs
from cogent.core.alignment import DenseAlignment
from qiime.filter_alignment import apply_lane_mask, apply_gap_filter,\
 apply_lane_make_and_gap_filters

class FilterAlignmentTests(TestCase):

    def setUp(self):
        """Init variables for the tests """
        self.aln1 = LoadSeqs(data=[\
         ('s1','ACC--T'),\
         ('s2','AC---T'),\
         ('s3','TCT--T'),\
         ('s4','ACG--T'),\
         ('s5','---A--'),\
         ],aligned=DenseAlignment)
         
    def test_apply_lane_make_and_gap_filters(self):
        """apply_lane_make_and_gap_filters: functions as expected
        """
        lm = '111111'
        expected = self.aln1
        self.assertEqual(apply_lane_make_and_gap_filters(\
         self.aln1,lm,1.0),expected)
         
        lm = None
        expected = self.aln1
        self.assertEqual(apply_lane_make_and_gap_filters(\
         self.aln1,lm,1.0),expected)
         
        # gap filter only
        lm = '111111'
        expected = LoadSeqs(data=[\
         ('s1','ACC-T'),\
         ('s2','AC--T'),\
         ('s3','TCT-T'),\
         ('s4','ACG-T'),\
         ('s5','---A-'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_lane_make_and_gap_filters(\
         self.aln1,lm),expected)
         
        # lm filter only
        lm = '011111'
        expected = LoadSeqs(data=[\
         ('s1','CC--T'),\
         ('s2','C---T'),\
         ('s3','CT--T'),\
         ('s4','CG--T'),\
         ('s5','--A--'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_lane_make_and_gap_filters(\
         self.aln1,lm,1.0),expected)
         
        # gap and lm filter only
        lm = '011111'
        expected = LoadSeqs(data=[\
         ('s1','CC-T'),\
         ('s2','C--T'),\
         ('s3','CT-T'),\
         ('s4','CG-T'),\
         ('s5','--A-'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_lane_make_and_gap_filters(\
         self.aln1,lm),expected)

    def test_apply_lane_mask(self):
        """ apply_lane_mask: functions as expected with varied lane masks
        """
        lm1 = '111111'
        expected = self.aln1
        self.assertEqual(apply_lane_mask(self.aln1,lm1),expected)        
        
        lm2 = '000000'
        expected = LoadSeqs(data=[\
         ('s1',''),\
         ('s2',''),\
         ('s3',''),\
         ('s4',''),\
         ('s5',''),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_lane_mask(self.aln1,lm2),expected)
        
        lm3 = '101010'
        expected = LoadSeqs(data=[\
         ('s1','AC-'),\
         ('s2','A--'),\
         ('s3','TT-'),\
         ('s4','AG-'),\
         ('s5','---'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_lane_mask(self.aln1,lm3),expected)
        
        lm4 = '000111'
        expected = LoadSeqs(data=[\
         ('s1','--T'),\
         ('s2','--T'),\
         ('s3','--T'),\
         ('s4','--T'),\
         ('s5','A--'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_lane_mask(self.aln1,lm4),expected)
        
    def test_apply_gap_filter(self):
        """ apply_gap_filter: functions as expected with varied allowed_gap_frac
        """
        expected = self.aln1
        self.assertEqual(apply_gap_filter(self.aln1,1.0),expected)
        
        expected = LoadSeqs(data=[\
         ('s1','ACC-T'),\
         ('s2','AC--T'),\
         ('s3','TCT-T'),\
         ('s4','ACG-T'),\
         ('s5','---A-'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_gap_filter(self.aln1),expected)
        
        expected = LoadSeqs(data=[\
         ('s1','ACCT'),\
         ('s2','AC-T'),\
         ('s3','TCTT'),\
         ('s4','ACGT'),\
         ('s5','----'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_gap_filter(self.aln1,0.75),expected)
           
        expected = LoadSeqs(data=[\
         ('s1','ACCT'),\
         ('s2','AC-T'),\
         ('s3','TCTT'),\
         ('s4','ACGT'),\
         ('s5','----'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_gap_filter(self.aln1,0.40),expected)
         
        expected = LoadSeqs(data=[\
         ('s1','ACT'),\
         ('s2','ACT'),\
         ('s3','TCT'),\
         ('s4','ACT'),\
         ('s5','---'),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_gap_filter(self.aln1,0.30),expected)
        
        expected = LoadSeqs(data=[\
         ('s1',''),\
         ('s2',''),\
         ('s3',''),\
         ('s4',''),\
         ('s5',''),\
         ],aligned=DenseAlignment)
        self.assertEqual(apply_gap_filter(self.aln1,0.10),expected)
        
    def test_apply_lane_make_and_gap_filters_alternate_alignment(self):
        """apply_lane_make_and_gap_filters: functions as expected with alt aln
        """
        aln = LoadSeqs(data=[\
         ('ACT009','AACT-'),\
         ('ACT019','AACT-'),\
         ('ACT011','-TCT-')],aligned=DenseAlignment)
        self.assertEqual(apply_lane_make_and_gap_filters(aln,None,1.0),aln)
        
        lm = '00111'
        expected = LoadSeqs(data=[\
         ('ACT009','CT'),\
         ('ACT019','CT'),\
         ('ACT011','CT')],aligned=DenseAlignment)
        self.assertEqual(apply_lane_make_and_gap_filters(aln,lm),expected)

if __name__ == "__main__":
    main()