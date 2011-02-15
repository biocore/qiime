#!/usr/bin/env python
# File created on 15 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.sort import (sort_sample_ids_by_mapping_value,
                       sort_fasta_by_abundance, natsort)

class SortTests(TestCase):
    
    def setUp(self):
        self.mapping_f1 = mapping_f1.split('\n')
        self.dirs_to_remove = []
        self.files_to_remove = []

    def tearDown(self):
        for dir in  self.dirs_to_remove:
            if exists(dir):
                rmdir(dir)
        remove_files(self.files_to_remove)
    
    def test_sort_sample_ids_by_mapping_value(self):
        """ sort_sample_ids_by_mapping_value functions as expected """
        actual = sort_sample_ids_by_mapping_value(mapping_file=self.mapping_f1,
                                         field='days_since_epoch',
                                         field_type_f=float)
        expected = zip(['NotInOtuTable','1','Z2','Z1','A'],
                       [0.0,5.7,10,23,400000])
        self.assertEqual(actual,expected)
        
    def test_sort_sample_ids_by_mapping_value_error(self):
        """ sort_sample_ids_by_mapping_value handles errors """
        self.assertRaises(ValueError,
                          sort_sample_ids_by_mapping_value,
                          mapping_file=self.mapping_f1,
                          field='years_since_spoch',
                          field_type_f=float)
                          
        self.assertRaises(ValueError,
                          sort_sample_ids_by_mapping_value,
                          mapping_file=self.mapping_f1,
                          field='Something',
                          field_type_f=float)
                          
    def test_sort_fasta_by_abundance(self):
      """sort_fasta_by_abundance functions as expected"""
      class FakeOutF(object):
          def __init__(self):
              self.s = ""
          def write(self,s):
              self.s += s

      actual = FakeOutF()
      expected = ""
      sort_fasta_by_abundance([],actual)
      self.assertEqual(actual.s,expected)

      # no sorting necessary
      actual = FakeOutF()
      expected1 = "\n".join(['>s1','ACCGT','>s2 comment','ATTGC',''])
      expected2 = "\n".join(['>s2 comment','ATTGC','>s1','ACCGT',''])
      sort_fasta_by_abundance(['>s1','ACCGT','>s2 comment','ATTGC'],actual)
      # order is unimportant here
      self.assertTrue(actual.s == expected1 or actual.s == expected2)

      # sorting necessary
      actual = FakeOutF()
      inseqs = ['>s1','ACCGT',
                 '>s2 comment','ATTGC',
                 '>s3 blah','ATTGC']
      expected = "\n".join(['>s2 comment','ATTGC',
                            '>s3 blah','ATTGC',
                            '>s1','ACCGT',''])
      sort_fasta_by_abundance(inseqs,actual)
      self.assertEqual(actual.s,expected)
      

    def test_natsort(self):
        """natsort should perform numeric comparisons on strings"""
        # string with alpha and numerics sort correctly
        s = 'sample1 sample2 sample11 sample12'.split()
        self.assertEqual(natsort(s), 
            'sample1 sample2 sample11 sample12'.split())
        s.reverse()
        self.assertEqual(natsort(s), 
            'sample1 sample2 sample11 sample12'.split())
        self.assertEqual(natsort(list('cba321')),list('123abc'))
        
        # strings with alpha only sort correctly
        self.assertEqual(natsort(list('cdba')),list('abcd'))
        
        # string of ints sort correctly
        self.assertEqual(natsort(['11','2','1','0']),
                                 ['0','1','2','11'])
                                 
        # strings of floats sort correctly
        self.assertEqual(natsort(['1.11','1.12','1.00','0.009']),
                                 ['0.009','1.00','1.11','1.12'])
          
          

mapping_f1 = """#SampleID\tSomething\tdays_since_epoch
Z1\t42\t23
Z2\thello\t10
A\t4\t400000
1\tr\t5.7
NotInOtuTable\tf\t0"""

if __name__ == "__main__":
    main()