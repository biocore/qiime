#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from convert_SampleMapping_to_OTUtable import SampleMapping_To_OTUtable, \
    parse_sample_mapping

class TestCode(TestCase):
    """Tests of the methods in convert_SampleMapping_to_OTUtable"""

    def setUp(self):
        """sets up variables for testing"""

        self.SampleMapping = """OTU1\tsample1\t3
OTU1\tsample3\t2
OTU2\tsample1\t1
OTU2\tsample2\t2""".split('\n')

    def test_parse_sample_mapping(self):
        """parse_sample_mapping works"""
        lines = self.SampleMapping
        OTU_sample_info, all_sample_names = parse_sample_mapping(lines)
        self.assertEqual(OTU_sample_info, {'OTU2': {'sample1': '1', 'sample3': '0', 'sample2': '2'}, 'OTU1': {'sample1': '3', 'sample3': '2', 'sample2': '0'}})
        
        self.assertEqual(all_sample_names, set(['sample1', 'sample3', 'sample2']))

    def test_SampleMapping_To_OTUtable(self):
        """SampleMapping_To_OTUtable works"""
        lines = self.SampleMapping
        result = SampleMapping_To_OTUtable(lines)
        self.assertEqual(result, ['#Full OTU Counts', '#OTU ID\tsample1\tsample2\tsample3', 'OTU2\t1\t2\t0', 'OTU1\t3\t0\t2'])

if __name__ == '__main__':
    main()
        
