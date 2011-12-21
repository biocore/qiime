#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from qiime.parse import parse_distmat
from qiime.filter import (filter_fasta,
                          filter_samples_from_distance_matrix,
                          negate_tips_to_keep)

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
        self.input_dm1 = input_dm1.split('\n')
        self.expected_dm1a = expected_dm1a.split('\n')
        self.expected_dm1b = expected_dm1b.split('\n')
        
        
    def tearDown(self):
        pass
        
    def test_negate_tips_to_keep(self):
        """ negate_tips_to_keep functions as expected """
        t = DndParser("((S5:0.00014,S7:0.00015)0.752:0.45762,(S3:0.00014,"
         "seq6:0.00014)0.180:0.00015,(Seq1:0.00014,s2:0.00014)0.528:1.0466);")
        
        tips_to_keep = ["S5","Seq1","s2"]
        expected = ["S7","S3","seq6"]
        self.assertEqualItems(negate_tips_to_keep(tips_to_keep,t),expected)
        
        tips_to_keep = ["S5","Seq1"]
        expected = ["S7","S3","seq6","s2"]
        self.assertEqualItems(negate_tips_to_keep(tips_to_keep,t),expected)
        
        tips_to_keep = []
        expected = ["S7","S3","seq6","s2","S5","Seq1"]
        self.assertEqualItems(negate_tips_to_keep(tips_to_keep,t),expected)
        
        tips_to_keep = ["S7","S3","seq6","s2","S5","Seq1"]
        expected = []
        self.assertEqualItems(negate_tips_to_keep(tips_to_keep,t),expected)
        
        
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
    
    def test_filter_samples_from_distance_matrix(self):
        """filter_samples_from_distance_matrix functions as expected """
        actual = filter_samples_from_distance_matrix(parse_distmat(self.input_dm1),
                                               ["GHI blah","XYZ"])
        self.assertEqual(actual,expected_dm1a)
        actual = filter_samples_from_distance_matrix(parse_distmat(self.input_dm1),
                                              ["GHI","DEF"])
        self.assertEqual(actual,expected_dm1b)
        
        
    def test_filter_samples_from_distance_matrix_file_input(self):
        """filter_samples_from_distance_matrix handles file input """
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                               ["GHI blah","XYZ"])
        self.assertEqual(actual,expected_dm1a)
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                              ["GHI","DEF"])
        self.assertEqual(actual,expected_dm1b)

    def test_filter_samples_from_distance_matrix_negate(self):
        """filter_samples_from_distance_matrix functions w negate """
        actual = filter_samples_from_distance_matrix(
          parse_distmat(self.input_dm1),
          ["ABC blah","DEF"],
          negate=True)
        self.assertEqual(actual,expected_dm1a)
        actual = filter_samples_from_distance_matrix(\
         parse_distmat(self.input_dm1),
         ["ABC","XYZ"],
         negate=True)
        self.assertEqual(actual,expected_dm1b)
        

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

input_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes
z\t0\t1\t0\t1\tNone
""" % __version__

expected_otu_table1a = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes
z\t0\t1\t0\t1\tNone""" % __version__

expected_otu_table1b = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
1\t1\t0\t0\t0\tNone
x\t0\t0\t3\t0\tBacteria;Bacteroidetes""" % __version__

input_seqs_to_discard1 = """x
1 some comment
42 not a real otu id"""

expected_otu_table1c = """# QIIME v%s OTU table
#OTU ID\tABC\tXYZ\tConsensus Lineage
0\t1\t0\tBacteria;Firmicutes
1\t1\t0\tNone
z\t0\t1\tNone""" % __version__

expected_otu_table1d = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tConsensus Lineage
0\t1\t1\t0\tBacteria;Firmicutes
1\t1\t0\t0\tNone
x\t0\t0\t3\tBacteria;Bacteroidetes
z\t0\t1\t0\tNone""" % __version__

input_dm1 = """\tABC\tDEF\tGHI\tXYZ
ABC\t0.0\t0.75\t0.00\t0.0063
DEF\t0.75\t0.0\t0.01\t0.65
GHI\t0.00\t0.01\t0.0\t1.0
XYZ\t0.0063\t0.065\t1.0\t0.0"""

expected_dm1a = """\tABC\tDEF
ABC\t0.0\t0.75
DEF\t0.75\t0.0"""

expected_dm1b = """\tABC\tXYZ
ABC\t0.0\t0.0063
XYZ\t0.0063\t0.0"""

input_otu_table2 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes;Something;Something else
1\t1\t0\t0\t0\tNone
0\t1\t1\t0\t2\tBacteria;Firmicutes;Something;Nothing
x\t0\t0\t3\t0\tBacteria;Bacteroidetes
z\t0\t1\t0\t1\tNone;Firmicutes
z\t0\t1\t0\t1\tArchaea""" % __version__

bacteria_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes;Something;Something else
0\t1\t1\t0\t2\tBacteria;Firmicutes;Something;Nothing
x\t0\t0\t3\t0\tBacteria;Bacteroidetes""" % __version__

archaea_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
z\t0\t1\t0\t1\tArchaea""" % __version__

none_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
1\t1\t0\t0\t0\tNone
z\t0\t1\t0\t1\tNone;Firmicutes""" % __version__

firmicutes_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t0\tBacteria;Firmicutes;Something;Something else
0\t1\t1\t0\t2\tBacteria;Firmicutes;Something;Nothing""" % __version__

bacteroidetes_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
x\t0\t0\t3\t0\tBacteria;Bacteroidetes""" % __version__

nothing_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
0\t1\t1\t0\t2\tBacteria;Firmicutes;Nothing""" % __version__

none_no_l2_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
1\t1\t0\t0\t0\tNone""" % __version__

none_firmicutes_otu_table1 = """# QIIME v%s OTU table
#OTU ID\tABC\tDEF\tGHI\tXYZ\tConsensus Lineage
z\t0\t1\t0\t1\tNone;Firmicutes""" % __version__

if __name__ == "__main__":
    main()