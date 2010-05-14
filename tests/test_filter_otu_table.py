#!/usr/bin/env python
#file test_filter_otu_table.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Release"

from sys import argv
from string import strip
from cogent.util.unit_test import TestCase, main
from numpy import array
from qiime.filter_otu_table import (strip_quotes,split_tax,
                _filter_table_samples)

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.otu_fname='otu_table.txt'
        self.min_count=1
        self.min_samples=2
        self.include_taxonomy='Root;Bacteria'
        self.exclude_taxonomy='Root;Archaea'
        self.dir_path='./'
        self.taxon1='"Root;Bacteria"'
        self.taxon2='Root;Archaea'

    def test_strip_quotes(self):
        """strip_quotes: removes leading and trailing quotes from tax string"""

        self.sample_to_extract='SampleA,SampleB'
        exp1='Root;Bacteria'
        exp2='Root;Archaea'
        
        obs1=strip_quotes(self.taxon1)
        obs2=strip_quotes(self.taxon2)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)

    def test_split_tax(self):
        """split_tax: Splits the tax string on comma and semicolon"""

        exp1=['"Root','Bacteria"']
        exp2=['Root','Archaea']
        
        obs1=split_tax(self.taxon1)
        obs2=split_tax(self.taxon2)
        
        self.assertEqual(obs1,exp1)
        self.assertEqual(obs2,exp2)
       
    def test_filter_table_samples(self):
        """_filter_table_samples removes small samples from OTU table
        """

        otu_table = """#Full OTU Counts
#OTU ID\tsample1\tsample2\tsample3
0\t0\t2\t0
1\t1\t0\t0
2\t1\t1\t1""".split('\n')
        result = _filter_table_samples(otu_table, 2)
        self.assertEqual(result, "#Full OTU Counts\n#OTU ID\tsample1\tsample2\n0\t0\t2\n1\t1\t0\n2\t1\t1")
        result = _filter_table_samples(otu_table, 1)
        self.assertEqual(result, '\n'.join(otu_table))
        result = _filter_table_samples(otu_table, 3)
        self.assertEqual(result, "#Full OTU Counts\n#OTU ID\tsample2\n0\t2\n1\t0\n2\t1")

#run tests if called from command line
if __name__ == "__main__":
    main()
