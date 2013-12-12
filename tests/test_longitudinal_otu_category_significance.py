#!/usr/bin/env python

"""Tests of code for perparing otu table for otu_category_significance of longitudinal studies
"""

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Catherine Lozupone", "Dan Knights", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.longitudinal_otu_category_significance import get_sample_individual_info, make_new_otu_counts, longitudinal_otu_table_conversion_wrapper
from numpy import array
from qiime.parse import parse_mapping_file
from biom.parse import parse_biom_table_str

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_get_sample_individual_info(self):
        """get_sample_individual_info works
        """
        mapping_lines = """#SampleID\tindividual\ttimepoint_zero\ttimepoint
AT0\tA\t1\t0
AT1\tA\t0\t1
AT2\tA\t0\t2
BT0\tB\t1\t0
BT1\tB\t0\t1
BT2\tB\t0\t2
""".split('\n')
        mapping_data, header, comments = parse_mapping_file(mapping_lines)
        samples_from_subject, samples_to_subtract = \
            get_sample_individual_info(mapping_data, header, 'individual', \
                'timepoint_zero')
        self.assertEqual(samples_from_subject, {'BT1': ['BT0', 'BT1', 'BT2'], 'BT0': ['BT0', 'BT1', 'BT2'], 'BT2': ['BT0', 'BT1', 'BT2'], 'AT2': ['AT0', 'AT1', 'AT2'], 'AT0': ['AT0', 'AT1', 'AT2'], 'AT1': ['AT0', 'AT1', 'AT2']})
        self.assertEqual(samples_to_subtract, {'BT1': 'BT0', 'BT0': 'BT0', 'BT2': 'BT0', 'AT2': 'AT0', 'AT0': 'AT0', 'AT1': 'AT0'})
        
    def test_make_new_otu_counts(self):
        """make_new_otu_counts works
        """
        mapping_lines = """#SampleID\tindividual\ttimepoint_zero\ttimepoint
AT0\tA\t1\t0
AT1\tA\t0\t1
AT2\tA\t0\t2
BT0\tB\t1\t0
BT1\tB\t0\t1
BT2\tB\t0\t2
""".split('\n')
        mapping_data, header, comments = parse_mapping_file(mapping_lines)
        samples_from_subject, sample_to_subtract = \
            get_sample_individual_info(mapping_data, header, 'individual', \
            'timepoint_zero')
        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 0.5], [0, 1, 0.29999999999999999], [0, 2, 99.0],
        [0, 3, 0.20000000000000001], [1, 2, 99.0], [1, 4, 0.40000000000000002],
        [1, 5, 0.5], [1, 6, 0.59999999999999998], [2, 0, 0.10000000000000001],
        [2, 1, 0.40000000000000002], [2, 2, 99.0], [2, 3, 0.69999999999999996],
        [2, 4, 0.5], [2, 5, 0.59999999999999998], [2, 6, 0.80000000000000004],
        [3, 1, 0.10000000000000001], [3, 2, 99.0],
        [3, 4, 0.40000000000000002]], "columns": [{"id": "AT0", "metadata":
        null}, {"id": "AT1", "metadata": null}, {"id": "S1", "metadata": null},
        {"id": "AT2", "metadata": null}, {"id": "BT0", "metadata": null},
        {"id": "BT1", "metadata": null}, {"id": "BT2", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 7], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T21:35:19.499263", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        converted_otu_table_object = make_new_otu_counts(otu_table,
                sample_to_subtract, samples_from_subject)
        expected_otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id":
        "3", "metadata": null}], "format":
        "Biological Observation Matrix 0.9.0-dev", "data": [[0, 1,
        -0.20000000000000001], [0, 2, -0.29999999999999999],
        [0, 3, 999999999.0], [0, 4, 999999999.0], [0, 5, 999999999.0], [1, 0,
        999999999.0], [1, 1, 999999999.0], [1, 2, 999999999.0], [1, 4,
        0.10000000000000001], [1, 5, 0.20000000000000001], [2, 1,
        0.29999999999999999], [2, 2, 0.59999999999999998], [2, 4,
        0.10000000000000001], [2, 5, 0.29999999999999999], [3, 1,
        0.10000000000000001], [3, 4, -0.40000000000000002], [3, 5,
        -0.40000000000000002]], "columns": [{"id": "AT0", "metadata": null},
        {"id": "AT1", "metadata": null}, {"id": "AT2", "metadata": null},
        {"id": "BT0", "metadata": null}, {"id": "BT1", "metadata": null},
        {"id": "BT2", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2570", "matrix_type": "sparse",
        "shape": [4, 6], "format_url":
        "http://biom-format.org",
        "date": "2011-12-21T21:43:06.809380", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        expected_otu_table_object = parse_biom_table_str(
                expected_otu_table_str)
        self.assertEqual(converted_otu_table_object.ObservationIds,
                         expected_otu_table_object.ObservationIds)
        self.assertEqual(converted_otu_table_object.SampleIds,
                         expected_otu_table_object.SampleIds)
        self.assertEqual(converted_otu_table_object.SampleMetadata,
                         expected_otu_table_object.SampleMetadata)
        self.assertEqual(converted_otu_table_object.ObservationMetadata,
                         expected_otu_table_object.ObservationMetadata)
        self.assertFloatEqual(sorted(converted_otu_table_object._data.items()),
                         sorted(expected_otu_table_object._data.items()))

    def test_longitudinal_otu_table_conversion_wrapper(self):
        """londitudinal_otu_table_conversion_wrapper works
        """
        mapping_lines = """#SampleID\tindividual\ttimepoint_zero\ttimepoint
AT0\tA\t1\t0
AT1\tA\t0\t1
AT2\tA\t0\t2
BT0\tB\t1\t0
BT1\tB\t0\t1
BT2\tB\t0\t2
""".split('\n')
        category_mapping = parse_mapping_file(mapping_lines)
        otu_table = """{"rows": [{"id": "0", "metadata": null}, {"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3", "metadata": null}, {"id": "4", "metadata": null}], "format": "Biological Observation Matrix 1.0.0", "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 3.0], [1, 3, 1.0], [1, 4, 2.0], [1, 5, 3.0], [2, 0, 1.0], [2, 1, 2.0], [2, 2, 3.0], [2, 4, 1.0], [2, 5, 2.0], [3, 0, 2.0], [3, 1, 4.0], [3, 2, 6.0], [3, 4, 1.0], [3, 5, 2.0], [4, 0, 3.0], [4, 1, 2.0], [4, 2, 1.0], [4, 3, 6.0], [4, 4, 4.0], [4, 5, 2.0]], "columns": [{"id": "AT0", "metadata": null}, {"id": "AT1", "metadata": null}, {"id": "AT2", "metadata": null}, {"id": "BT0", "metadata": null}, {"id": "BT1", "metadata": null}, {"id": "BT2", "metadata": null}], "generated_by": "BIOM-Format 1.0.0-dev", "matrix_type": "sparse", "shape": [5, 6], "format_url": "http://biom-format.org", "date": "2012-08-01T09:14:03.574451", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""
    
        otu_table = parse_biom_table_str(otu_table)
        new_otu_table = longitudinal_otu_table_conversion_wrapper(otu_table,\
            category_mapping, 'individual', 'timepoint_zero')
        new_otu_table = str(new_otu_table).split('\n')
        self.assertEqual(new_otu_table[0], "# Constructed from biom file")
        data_line1 = new_otu_table[2].split('\t')
        self.assertFloatEqual(float(data_line1[0]), 0.0) 
        #sets the reference to 0
        self.assertFloatEqual(float(data_line1[1]), 0.0)
        #subtracts values from same individual from the reference
        self.assertFloatEqual(float(data_line1[2]), 0.05714286)
        #sets to ignore number when not observed across a person
        self.assertFloatEqual(float(data_line1[4]), 999999999.0)


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
