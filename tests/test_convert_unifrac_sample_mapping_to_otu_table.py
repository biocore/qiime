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
    sample_mapping_to_otu_table, sample_mapping_to_biom_table)

class ConvertUnifracSampleMappingToOtuTableTests(TestCase):
    """"""

    def setUp(self):
        """"""
        self.SampleMapping = ["OTU1\tsample1\t3", "OTU1\tsample3\t2", \
        "OTU2\tsample1\t1", "OTU2\tsample2\t2"]

        self.SampleMapping2 = ["OTU1\tsample1", "OTU1\tsample3", \
        "OTU2\tsample1", "OTU2\tsample2"]

    def test_sample_mapping_to_otu_table(self):
        """sample_mapping_to_otu_table works"""
        lines = self.SampleMapping
        result = sample_mapping_to_otu_table(lines)
        self.assertEqual(result, ['#Full OTU Counts',\
         '#OTU ID\tsample1\tsample2\tsample3', 'OTU2\t1\t2\t0', \
        'OTU1\t3\t0\t2'])

    def test_sample_mapping_to_biom_table(self):
        """sample_mapping_to_biom_table works"""
        lines = self.SampleMapping
        actual = sample_mapping_to_biom_table(lines)
        exp = table_factory(array([[3.,0.,2.],[1.,2.,0.]]),
                            ['sample1','sample2','sample3'],
                            ['OTU1','OTU2'])
        self.assertEqual(actual.sortBySampleId(), exp.sortBySampleId())

if __name__ =='__main__':
    main()