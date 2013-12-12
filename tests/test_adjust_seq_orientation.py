#!/usr/bin/env python
# File created on 07 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from unittest import TestCase, main
from cogent.parse.fasta import MinimalFastaParser
from qiime.adjust_seq_orientation import rc_fasta_lines, null_seq_desc_mapper,\
    append_rc
    
class AdjustSeqOrientationTests(TestCase):
    """ """
    
    def setUp(self):
        """ """
        self.fasta_lines1 = fasta_lines1.split('\n')
        self.fasta_lines1_mixed_case = fasta_lines1_mixed_case.split('\n')
        
        self.fasta_lines1_exp = list(MinimalFastaParser(
            fasta_lines1_exp.split('\n')))
        self.fasta_lines1_mixed_case_exp = list(MinimalFastaParser(
         fasta_lines1_mixed_case_exp.split('\n')))
        self.fasta_lines1_exp_null_desc_mapper = list(MinimalFastaParser(
            fasta_lines1_exp_null_desc_mapper.split('\n')))
            
    def test_rc_fasta_lines(self):
        """rc_fasta_lines: functions as expected w/ seq_id mapping 
        """
        self.assertEqual(list(rc_fasta_lines(self.fasta_lines1,append_rc)),
                         self.fasta_lines1_exp)
            
    def test_rc_fasta_lines_mixed_case(self):
        """rc_fasta_lines: functions with mixed cases in sequences
        """
        self.assertEqual(list(
         rc_fasta_lines(self.fasta_lines1_mixed_case,append_rc)),
         self.fasta_lines1_mixed_case_exp)
            
    def test_rc_fasta_lines_leave_seq_desc(self):
        """rc_fasta_lines: functions as expected w/o seq_id mapping 
        """
        self.assertEqual(list(
         rc_fasta_lines(self.fasta_lines1,null_seq_desc_mapper)),
                         self.fasta_lines1_exp_null_desc_mapper)


fasta_lines1 = """>s1 some description
AAATGGCGCGCG
>s2
TTATATCCGC
"""

fasta_lines1_mixed_case = """>s1 some description
aaatGGcgcgcg
>s2
ttatatccgc
"""

fasta_lines1_exp = """>s1 some description RC
CGCGCGCCATTT
>s2 RC
GCGGATATAA
"""

fasta_lines1_mixed_case_exp = """>s1 some description RC
CGCGCGCCATTT
>s2 RC
GCGGATATAA
"""

fasta_lines1_exp_null_desc_mapper = """>s1 some description
CGCGCGCCATTT
>s2
GCGGATATAA
"""



if __name__ == "__main__":
    main()