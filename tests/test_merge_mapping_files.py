#!/usr/bin/env python
# File created on 30 Nov 2009.
from __future__ import division
from qiime.merge_mapping_files import merge_mapping_files
from unittest import TestCase, main

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


class MergeMappingFilesTests(TestCase):
    """ Tests of the MergeMappingFiles script"""
    
    def setUp(self):
        """ """
        self.m1 = m1.split('\n')
        self.m2 = m2.split('\n')
        self.m3 = m3.split('\n')
        self.m1_m2_m3_exp = m1_m2_m3_exp.split('\n')
    
    def test_merge_mapping_file(self):
        """merge_mapping_file: functions with default parameters """
        actual = merge_mapping_files([self.m1,self.m2,self.m3])
        expected = self.m1_m2_m3_exp
        
        self.assertTrue(actual[0].startswith('#SampleID\tBarcodeSequence'))
        self.assertTrue(actual[0].endswith('Description'))
        
        actual.sort()
        expected.sort()
        
        for a,e in zip(actual,expected):
            a_fields = a.split('\t')
            e_fields = e.split('\t')
            a_fields.sort()
            e_fields.sort()
            self.assertEqual(a_fields,e_fields)
            
    def test_merge_mapping_file_dups(self):
        """merge_mapping_file: duplicate sample ids stay duplicate """
        actual = merge_mapping_files([self.m1,self.m2,self.m1])
        # length is number of samples plus header line
        self.assertEqual(len(actual),6)
            
    def test_merge_mapping_correct_num_records(self):
        """merge_mapping_file: num recs returned is correct with varied input 
        """
        # length is always number of samples plus 1 (for header line)
        actual = merge_mapping_files([self.m1])
        self.assertEqual(len(actual),3)
        actual = merge_mapping_files([self.m2])
        self.assertEqual(len(actual),2)
        actual = merge_mapping_files([self.m3])
        self.assertEqual(len(actual),4)
        
        actual = merge_mapping_files([self.m1,self.m2])
        self.assertEqual(len(actual),4)
        actual = merge_mapping_files([self.m1,self.m3])
        self.assertEqual(len(actual),6)
        actual = merge_mapping_files([self.m2,self.m3])
        self.assertEqual(len(actual),5)
        
        actual = merge_mapping_files([self.m1,self.m2,self.m3])
        self.assertEqual(len(actual),7)
        
        actual = merge_mapping_files([self.m3,self.m3,self.m3,self.m3])
        self.assertEqual(len(actual),13)
        


m1 = """#SampleID\tBarcodeSequence\tdata1\tdata2\tDescription
samp1_1\tAAAA\t42\t36.9\tsomething
samp1_2\tAAAA\t99\t22.5\t"nothing interesting"
"""

m2 = """#SampleID\tBarcodeSequence\tdata 3\tdata4\tdata5\tdata2\tDescription
samp2_1\tAAAA\tgreen\tsoil\t99.8\t44.5\tother
"""

m3 = """#SampleID\tBarcodeSequence\tdata1\tDescription
samp3_1\tAAAT\t8\tmisc1
samp3_2\tAAAG\t6\tmisc2
samp3_3\tAAAC\t7\tmisc3
"""

m1_m2_m3_exp = """#SampleID\tBarcodeSequence\tdata1\tdata2\tdata 3\tdata4\tdata5\tDescription
samp1_1\tAAAA\t42\t36.9\tno_data\tno_data\tno_data\tsomething
samp1_2\tAAAA\t99\t22.5\tno_data\tno_data\tno_data\t"nothing interesting"
samp2_1\tAAAA\tno_data\t44.5\tgreen\tsoil\t99.8\tother
samp3_1\tAAAT\t8\tno_data\tno_data\tno_data\tno_data\tmisc1
samp3_2\tAAAG\t6\tno_data\tno_data\tno_data\tno_data\tmisc2
samp3_3\tAAAC\t7\tno_data\tno_data\tno_data\tno_data\tmisc3"""

if __name__ == "__main__":
    main()