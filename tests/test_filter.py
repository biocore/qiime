#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from qiime.filter import filter_fasta

class fake_output_f():
    
    def __init__(self):
        self.s = ""
    
    def write(self,s):
        self.s += s
    
    def close(self):
        pass

class FilterTests(TestCase):
    
    def setUp(self):
        self.filter_fasta_expected1 = filter_fasta_expected1
        self.filter_fasta_expected2 = filter_fasta_expected2
        
    def tearDown(self):
        pass
        
    def test_filter_fasta(self):
        """filter_fasta functions as expected """
        input_seqs = [('Seq1 some comment','ACCTTGG'),
                      ('s2 some other comment','TTGG'),
                      ('S3','AAGGCCGG'),
                      ('S5 some comment','CGT'),
                      ('seq6 some other comment','AA'),
                      ('S7','T')]
        seqs_to_keep = {}.fromkeys(['Seq1',
                                    's2 some other comment',
                                    'S3 no comment'])

        actual = fake_output_f()
        filter_fasta(input_seqs,
                     actual,
                     seqs_to_keep,
                     negate=False)
        self.assertEqual(actual.s,self.filter_fasta_expected1)
        
        actual = fake_output_f()
        filter_fasta(input_seqs,
                     actual,
                     seqs_to_keep,
                     negate=True)
        self.assertEqual(actual.s,self.filter_fasta_expected2)

filter_fasta_expected1 = """>Seq1 some comment
ACCTTGG
>s2 some other comment
TTGG
>S3
AAGGCCGG
"""
filter_fasta_expected2 = """>S5 some comment
CGT
>seq6 some other comment
AA
>S7
T
"""

if __name__ == "__main__":
    main()