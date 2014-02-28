#!/usr/bin/env python
# File created on 30 Nov 2009.
from __future__ import division

from unittest import TestCase, main
from StringIO import StringIO

from qiime.util import MetadataMap
from qiime.merge_mapping_files import merge_mapping_files

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jesse Stombaugh", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

class MergeMappingFilesTests(TestCase):

    """ Tests of the MergeMappingFiles script"""
    def setUp(self):
        self.m1 = StringIO(m1)
        self.m2 = StringIO(m2)
        self.m3 = StringIO(m3)

        self.m1_dup_bad = StringIO(m1_dup_bad)
        self.m1_dup_good = StringIO(m1_dup_good)
    
    def test_merge_mapping_file(self):
        """merge_mapping_file: functions with default parameters """
        observed = merge_mapping_files([self.m1,self.m3])
        self.assertEqual(observed.__dict__, m1_m3_exp)
            
    def test_merge_mapping_file_different_no_data_value(self):
        """merge_mapping_file: functions with different no_data_value"""
        observed = merge_mapping_files([self.m1,self.m2], "TESTING_NA")
        self.assertEqual(observed.__dict__, m1_m2_exp)
            
    def test_merge_mapping_file_three_mapping_files(self):
        """merge_mapping_file: 3 mapping files"""
        observed = merge_mapping_files([self.m1,self.m2,self.m3])
        self.assertEqual(observed.__dict__, m1_m2_m3_exp)
            
    def test_merge_mapping_file_bad_duplicates(self):
        """merge_mapping_file: error raised when merging mapping files where same sample ids has different values """
        self.assertRaises(
            ValueError,
            merge_mapping_files,
            [self.m1,
             self.m1_dup_bad])

    def test_merge_mapping_file_good_duplicates(self):
        """merge_mapping_file: same sample ids merged correctly when they have mergable data
        """
        actual = merge_mapping_files([self.m1,self.m1_dup_good])

        self.assertEqual(actual.__dict__, m1_m1_dup_good_exp)

m1="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional1\tDescription
111111111\tAAAAAAAAAAAAAAA\tTTTTTTTTTTTTTTTTTTTT\tfirst1111\tTHE FIRST"""

# when merged with m1, this should raise a ValueError
m1_dup_bad="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional1\tDescription
111111111\tAAAAAAAAAAAAAAA\tTTTTTTTTTTTTTTTTTTTT\tdupe11111\tDUPE BAD"""

# when merged with m1, this should NOT raise a ValueError
m1_dup_good="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional9\tDescription
111111111\tAAAAAAAAAAAAAAA\tTTTTTTTTTTTTTTTTTTTT\tsomthing\tTHE FIRST"""

m2="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional2\tDescription
222222222\tGGGGGGGGGGGGGGG\tCCCCCCCCCCCCCCCCCCCC\tsecond222\tTHE SECOND"""

m3="""#SampleID\tBarcodeSequence\tLinkerPrimerSequence\toptional1\tDescription
333333333\tFFFFFFFFFFFFFFF\tIIIIIIIIIIIIIIIIIIII\tthird3333\tTHE THIRD"""

# the dict representation of the object with no_data_value set to TESTING_NA
m1_m2_exp={'_metadata': {'111111111': {'optional1': 'first1111', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'BarcodeSequence': 'AAAAAAAAAAAAAAA', 'Description': 'THE FIRST', 'optional2': None}, '222222222': {'optional1': None, 'LinkerPrimerSequence': 'CCCCCCCCCCCCCCCCCCCC', 'BarcodeSequence': 'GGGGGGGGGGGGGGG', 'Description': 'THE SECOND', 'optional2': 'second222'}}, 'no_data_value': 'TESTING_NA', 'Comments': []}

# dict representation, with no_data_value left as default ("no_data")
m1_m3_exp={'_metadata': {'111111111': {'optional1': 'first1111', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'BarcodeSequence': 'AAAAAAAAAAAAAAA', 'Description': 'THE FIRST'}, '333333333': {'optional1': 'third3333', 'LinkerPrimerSequence': 'IIIIIIIIIIIIIIIIIIII', 'BarcodeSequence': 'FFFFFFFFFFFFFFF', 'Description': 'THE THIRD'}}, 'no_data_value': 'no_data', 'Comments': []}

# dict representation, with no_data_value left as default ("no_data")
m1_m2_m3_exp={'_metadata': {'111111111': {'BarcodeSequence': 'AAAAAAAAAAAAAAA', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'optional1': 'first1111', 'Description': 'THE FIRST', 'optional2': None}, '333333333': {'BarcodeSequence': 'FFFFFFFFFFFFFFF', 'LinkerPrimerSequence': 'IIIIIIIIIIIIIIIIIIII', 'optional1': 'third3333', 'Description': 'THE THIRD', 'optional2': None}, '222222222': {'BarcodeSequence': 'GGGGGGGGGGGGGGG', 'LinkerPrimerSequence': 'CCCCCCCCCCCCCCCCCCCC', 'optional1': None, 'Description': 'THE SECOND', 'optional2': 'second222'}}, 'no_data_value': 'no_data', 'Comments': []}

# dict representation, m1 and m1_dup_bad should be mergeable
m1_m1_dup_good_exp = {'_metadata': {'111111111': {'optional9': 'somthing', 'LinkerPrimerSequence': 'TTTTTTTTTTTTTTTTTTTT', 'optional1': 'first1111', 'Description': 'THE FIRST', 'BarcodeSequence': 'AAAAAAAAAAAAAAA'}}, 'no_data_value': 'no_data', 'Comments': []}

if __name__ == "__main__":
    main()
