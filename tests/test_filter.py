#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_otu_table
from qiime.filter import (filter_fasta, 
                          filter_otus_from_otu_table,
                          filter_samples_from_otu_table,
                          filter_otu_table_to_n_samples,
                          filter_samples_from_distance_matrix)

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
        self.expected_otu_table1c = expected_otu_table1c.split('\n')
        self.expected_otu_table1d = expected_otu_table1d.split('\n')
        self.input_dm1 = input_dm1.split('\n')
        self.expected_dm1a = expected_dm1a.split('\n')
        self.expected_dm1b = expected_dm1b.split('\n')
        
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
        
    def test_filter_samples_from_otu_table(self):
        """filter_samples_from_otu_table functions as expected """
        
        actual = filter_samples_from_otu_table(self.input_otu_table1,
                                              ["DEF","GHI tfasd"])
        self.assertEqual(actual,self.expected_otu_table1c)
        
        # order of otu table is retained regardless of samples_to_keep order
        actual = filter_samples_from_otu_table(self.input_otu_table1,
                                               ["XYZ"])
        self.assertEqual(actual,self.expected_otu_table1d)
        
    def test_filter_samples_from_otu_table_negate(self):
        """filter_samples_from_otu_table functions w negate """
        actual = filter_samples_from_otu_table(self.input_otu_table1,
                                               ["ABC blah","XYZ"],
                                               negate=True)
        self.assertEqual(actual,self.expected_otu_table1c)
    
    def test_filter_otu_table_to_n_samples(self):
        """filter_otu_table_to_n_samples functions as expected """
        
        def check_no_zero_count_otus(data,otu_ids):
            for d,id_ in zip(data,otu_ids):
                self.assertTrue(sum(d)>0,
                 "OTU (%s) with total count=0 not " % id_ +\
                 "filtered from OTU table.")
        
        actual = filter_otu_table_to_n_samples(self.input_otu_table1,1)
        sample_ids,otu_ids,data,taxa = parse_otu_table(actual)
        self.assertTrue(len(sample_ids),1)
        check_no_zero_count_otus(data,otu_ids)
        
        actual = filter_otu_table_to_n_samples(self.input_otu_table1,2)
        sample_ids,otu_ids,data,taxa = parse_otu_table(actual)
        self.assertTrue(len(sample_ids),2)
        check_no_zero_count_otus(data,otu_ids)
        
        actual = filter_otu_table_to_n_samples(self.input_otu_table1,3)
        sample_ids,otu_ids,data,taxa = parse_otu_table(actual)
        self.assertTrue(len(sample_ids),3)
        check_no_zero_count_otus(data,otu_ids)
        
        actual = filter_otu_table_to_n_samples(self.input_otu_table1,4)
        sample_ids,otu_ids,data,taxa = parse_otu_table(actual)
        self.assertTrue(len(sample_ids),4)
        self.assertEqualItems(sample_ids,["ABC","DEF","GHI","XYZ"])
        self.assertEqual(data.shape,(4,4))
        check_no_zero_count_otus(data,otu_ids)
        
        actual = filter_otu_table_to_n_samples(self.input_otu_table1,5)
        sample_ids,otu_ids,data,taxa = parse_otu_table(actual)
        self.assertTrue(len(sample_ids),4)
        self.assertEqualItems(sample_ids,["ABC","DEF","GHI","XYZ"])
        self.assertEqual(data.shape,(4,4))
        check_no_zero_count_otus(data,otu_ids)
        
        # different samples selected on different runs
        r = []
        for i in range(25):
            actual = filter_otu_table_to_n_samples(self.input_otu_table1,1)
            sample_ids,otu_ids,data,taxa = parse_otu_table(actual)
            r.append(tuple(sample_ids))
        self.assertTrue(len({}.fromkeys(r)) > 1)

    def test_filter_otu_table_to_n_samples_handles_bad_n(self):
        """filter_otu_table_to_n_samples handles bad n """
        # ValueError if n < 1
        actual = self.assertRaises(ValueError,
                                   filter_otu_table_to_n_samples,
                                   self.input_otu_table1,
                                   0)
        actual = self.assertRaises(ValueError,
                                   filter_otu_table_to_n_samples,
                                   self.input_otu_table1,
                                   -1)
        
    def test_filter_samples_from_distance_matrix(self):
        """filter_samples_from_distance_matrix functions as expected """
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                               ["GHI blah","XYZ"])
        self.assertEqual(actual,expected_dm1a)
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                              ["GHI","DEF"])
        self.assertEqual(actual,expected_dm1b)

    def test_filter_samples_from_distance_matrix_negate(self):
        """filter_samples_from_distance_matrix functions w negate """
        actual = filter_samples_from_distance_matrix(self.input_dm1,
                                               ["ABC blah","DEF"],
                                               negate=True)
        self.assertEqual(actual,expected_dm1a)
        actual = filter_samples_from_distance_matrix(self.input_dm1,
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

if __name__ == "__main__":
    main()