#!/usr/bin/env python
# File created on 08 Jun 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from unittest import TestCase, main
from biom.parse import parse_biom_table
from qiime.core_microbiome import (core_observations_across_sample_ids)


class ComputeCoreMicrobiomeTests(TestCase):

    """ """

    def setUp(self):
        """ """
        self.otu_table_data1 = parse_biom_table(otu_table1)
        self.otu_table_data2 = parse_biom_table(otu_table2)

    def test_core_observations_across_sample_ids(self):
        """ core_observations_across_sample_ids functions as expected
        """
        actual = core_observations_across_sample_ids(self.otu_table_data1,
                                                     ["S1", "s2"],
                                                     fraction_for_core=1.)
        expected = ['o1', 'o5']
        self.assertEqual(actual, expected)

        # fraction_for_core = 0.5
        actual = core_observations_across_sample_ids(self.otu_table_data1,
                                                     ["S1", "s2"],
                                                     fraction_for_core=0.5)
        expected = ['o1', 'o3', 'o5']
        self.assertEqual(actual, expected)

    def test_core_observations_across_sample_ids_invalid(self):
        """ core_observations_across_sample_ids handles invalid input as expected
        """
        self.assertRaises(ValueError,
                          core_observations_across_sample_ids,
                          self.otu_table_data1,
                          ["S1", "s2"],
                          fraction_for_core=1.001)

        self.assertRaises(ValueError,
                          core_observations_across_sample_ids,
                          self.otu_table_data1,
                          ["S1", "s2"],
                          fraction_for_core=-0.001)

    def test_core_observations_across_sample_ids_no_core(self):
        """core_observations_across_sample_ids handles filtering all obs
        """
        actual = core_observations_across_sample_ids(self.otu_table_data2,
                                                     ["S1", "s2", "s3", "s4"],
                                                     fraction_for_core=1.)
        expected = []
        self.assertEqual(actual, expected)


otu_table1 = """{"rows": [{"id": "o1", "metadata": {"OTUMetaData": "Eukarya;Human"}}, {"id": "o2", "metadata": {"OTUMetaData": "Eukarya;Moose"}}, {"id": "o3", "metadata": {"OTUMetaData": "Eukarya;Galapagos Tortoise"}}, {"id": "o4", "metadata": {"OTUMetaData": "Eukarya;Bigfoot"}}, {"id": "o5", "metadata": {"OTUMetaData": "Eukarya;Chicken"}}], "format": "Biological Observation Matrix 0.9.3", "data": [[0, 0, 105.0], [0, 1, 42.0], [0, 2, 99.0], [0, 3, 60000.0], [1, 2, 9.0], [1, 3, 99.0], [2, 0, 45.0], [4, 0, 1.0], [4, 1, 2.0], [4, 3, 3.0]], "columns": [{"id": "S1", "metadata": null}, {"id": "s2", "metadata": null}, {"id": "s3", "metadata": null}, {"id": "s4", "metadata": null}], "generated_by": "BIOM-Format 0.9.3", "matrix_type": "sparse", "shape": [5, 4], "format_url": "http://biom-format.org", "date": "2012-06-08T14:42:46.058411", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

otu_table2 = """{"rows": [{"id": "o1", "metadata": null}, {"id": "o2", "metadata": null}, {"id": "o3", "metadata": null}, {"id": "o4", "metadata": null}, {"id": "o5", "metadata": null}], "format": "Biological Observation Matrix 0.9.3", "data": [[0, 0, 105.0], [0, 1, 42.0], [0, 2, 99.0], [1, 2, 9.0], [1, 3, 99.0], [2, 0, 45.0], [4, 0, 1.0], [4, 1, 2.0], [4, 3, 3.0]], "columns": [{"id": "S1", "metadata": null}, {"id": "s2", "metadata": null}, {"id": "s3", "metadata": null}, {"id": "s4", "metadata": null}], "generated_by": "BIOM-Format 0.9.3", "matrix_type": "sparse", "shape": [5, 4], "format_url": "http://biom-format.org", "date": "2012-06-08T14:43:27.964500", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

if __name__ == "__main__":
    main()
