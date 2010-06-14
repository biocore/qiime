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
from qiime.filter import filter_fasta, filter_otus_from_otu_table

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
        self.input_otu_table1 = input_otu_table1.split('\n')
        self.input_seqs_to_discard1 = input_seqs_to_discard1.split('\n')
        self.expected_otu_table1a = expected_otu_table1a.split('\n')
        self.expected_otu_table1b = expected_otu_table1b.split('\n')
        
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
        
    def test_filter_otus_from_otu_table(self):
        """filter_otus_from_otu_table: functions as expected
        """
        actual = filter_otus_from_otu_table(self.input_otu_table1,
                                            self.input_seqs_to_discard1)
        self.assertEqual(actual,self.expected_otu_table1a)
        
    def test_filter_otus_from_otu_table_negate(self):
        """filter_otus_from_otu_table: functions as expected with negate=True
        """
        actual = filter_otus_from_otu_table(self.input_otu_table1,
                                            self.input_seqs_to_discard1,
                                            negate=True)
        self.assertEqual(actual,self.expected_otu_table1b)
        
        

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

input_otu_table1 = """#Full OTU Counts
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes
z\t0\t1\t0\t1\tNone
"""

expected_otu_table1a = """#Full OTU Counts
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
z\t0\t1\t0\t1\tNone"""

expected_otu_table1b = """#Full OTU Counts
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes"""

input_seqs_to_discard1 = """x
1 some comment
42 not a real otu id"""

if __name__ == "__main__":
    main()