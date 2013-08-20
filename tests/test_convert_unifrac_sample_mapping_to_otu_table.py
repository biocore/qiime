#! /usr/bin/env python

__author__ = "Cathy Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone", "Greg Caporaso",
                "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Cathy Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase,main
from numpy import array
from biom.table import table_factory
from qiime.convert_unifrac_sample_mapping_to_otu_table import (
    parse_sample_mapping, sample_mapping_to_otu_table,
    sample_mapping_to_biom_table)

class ConvertUnifracSampleMappingToOtuTableTests(TestCase):
    """"""

    def setUp(self):
        """"""
        self.SampleMapping = ["OTU1\tsample1\t3", "OTU1\tsample3\t2",
        "OTU2\tsample1\t1", "OTU2\tsample2\t2"]

        self.SampleMappingNoMIENS = ["OTU1\tsample_1\t3", "OTU1\tsample#3\t2",
        "OTU2\tsample_1\t1", "OTU2\tsample@2\t2"]

        self.SampleMapping2 = ["OTU1\tsample1", "OTU1\tsample3", \
        "OTU2\tsample1", "OTU2\tsample2"]

    def test_parse_sample_mapping(self):
        """parse_sample_mapping works"""
        lines = self.SampleMapping
        obs_OTU_sample_info, obs_all_sample_names = parse_sample_mapping(lines)
        exp_OTU_sample_info = {
            'OTU2': {'sample1': '1', 'sample3': '0', 'sample2': '2'},
            'OTU1': {'sample1': '3', 'sample3': '2', 'sample2': '0'}
        }
        exp_all_sample_names = set(['sample1', 'sample3', 'sample2'])
        self.assertEqual(obs_OTU_sample_info, exp_OTU_sample_info)
        self.assertEqual(obs_all_sample_names, exp_all_sample_names)
        #test that it corrects no MIENS compliant sample ids
        lines = self.SampleMappingNoMIENS
        obs_OTU_sample_info, obs_all_sample_names = parse_sample_mapping(lines)
        exp_OTU_sample_info = {
            'OTU2': {'sample.1': '1', 'sample.3': '0', 'sample.2': '2'},
            'OTU1': {'sample.1': '3', 'sample.3': '2', 'sample.2': '0'}
        }
        exp_all_sample_names = set(['sample.1', 'sample.3', 'sample.2'])
        self.assertEqual(obs_OTU_sample_info, exp_OTU_sample_info)
        self.assertEqual(obs_all_sample_names, exp_all_sample_names)
        #test that it works if no sample counts in there
        lines = self.SampleMapping2
        obs_OTU_sample_info, obs_all_sample_names = parse_sample_mapping(lines)
        exp_OTU_sample_info = {
            'OTU2': {'sample1': '1', 'sample3': '0', 'sample2': '1'},
            'OTU1': {'sample1': '1', 'sample3': '1', 'sample2': '0'}
        }
        exp_all_sample_names = set(['sample1', 'sample3', 'sample2'])
        self.assertEqual(obs_OTU_sample_info, exp_OTU_sample_info)
        self.assertEqual(obs_all_sample_names, exp_all_sample_names)

    def test_sample_mapping_to_otu_table(self):
        """sample_mapping_to_otu_table works"""
        lines = self.SampleMapping
        result = sample_mapping_to_otu_table(lines)
        self.assertEqual(result, ['#Full OTU Counts',\
         '#OTU ID\tsample1\tsample2\tsample3', 'OTU2\t1\t2\t0', \
        'OTU1\t3\t0\t2'])

        lines = self.SampleMappingNoMIENS
        result = sample_mapping_to_otu_table(lines)
        self.assertEqual(result, ['#Full OTU Counts',\
         '#OTU ID\tsample.1\tsample.2\tsample.3', 'OTU2\t1\t2\t0', \
        'OTU1\t3\t0\t2'])

    def test_sample_mapping_to_biom_table(self):
        """sample_mapping_to_biom_table works"""
        lines = self.SampleMapping
        actual = sample_mapping_to_biom_table(lines)
        exp = table_factory(array([[3.,0.,2.],[1.,2.,0.]]),
                            ['sample1','sample2','sample3'],
                            ['OTU1','OTU2'])
        self.assertEqual(actual.sortBySampleId(), exp.sortBySampleId())

        lines = self.SampleMappingNoMIENS
        actual = sample_mapping_to_biom_table(lines)
        exp = table_factory(array([[3.,0.,2.],[1.,2.,0.]]),
                            ['sample.1','sample.2','sample.3'],
                            ['OTU1','OTU2'])
        self.assertEqual(actual.sortBySampleId(), exp.sortBySampleId())


if __name__ =='__main__':
    main()