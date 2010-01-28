#!/usr/bin/env python
#file test_filter_otu_table.py

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"

from qiime.parse import parse_otus
from sys import argv
from string import strip
from cogent.util.unit_test import TestCase, main
from numpy import array
from qiime.filter_otu_table import (strip_quotes,split_tax,process_options)

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
        
    def test_process_options(self):
        """process_options: converts all options to a dict for passing"""
        
        exp={}
        exp['otu_file'] = 'otu_table.txt'
        exp['min_otu_count'] = 1
        exp['min_otu_samples'] = 2
        exp['included_taxa']=set(['Root','Bacteria'])
        exp['excluded_taxa']=set(['Root','Archaea'])
        exp['filtered_otu_file']='.//otu_table_filtered.txt'
        
        obs=process_options(self)
        
        self.assertEqual(obs,exp)

        


#run tests if called from command line
if __name__ == "__main__":
    main()
