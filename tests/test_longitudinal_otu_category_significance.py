#!/usr/bin/env python

"""Tests of code for perparing otu table for otu_category_significance of longitudinal studies
"""

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Catherine Lozupone", "Dan Knights"] 
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Release"

from cogent.util.unit_test import TestCase, main
from qiime.longitudinal_otu_category_significance import get_sample_individual_info, make_new_otu_counts
from numpy import array
from qiime.parse import parse_otu_table, parse_mapping_file

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
        otu_lines = """# QIIME v1.2.0-dev OTU table
#OTU ID\tAT0\tAT1\tS1\tAT2\tBT0\tBT1\tBT2
0\t0.5\t0.3\t99\t0.2\t0.0\t0.0\t0.0
1\t0.0\t0.0\t99\t0.0\t0.4\t0.5\t0.6
2\t0.1\t0.4\t99\t0.7\t0.5\t0.6\t0.8
3\t0.0\t0.1\t99\t0.0\t0.4\t0.0\t0.0
""".split('\n')
        otu_table = parse_otu_table(otu_lines, float)
        sample_ids, otu_ids, otu_counts, consensus = otu_table
        converted_otu_table = make_new_otu_counts(otu_ids, sample_ids, otu_counts, consensus, sample_to_subtract, samples_from_subject)
        converted_otu_table = converted_otu_table.split('\n')
        self.assertEqual(converted_otu_table[1], "#OTU ID\tAT0\tAT1\tAT2\tBT0\tBT1\tBT2")
        self.assertEqual(converted_otu_table[2], "0\t0.0\t-0.2\t-0.3\t999999999.0\t999999999.0\t999999999.0")
        self.assertEqual(converted_otu_table[3], "1\t999999999.0\t999999999.0\t999999999.0\t0.0\t0.1\t0.2")
        self.assertEqual(converted_otu_table[4], "2\t0.0\t0.3\t0.6\t0.0\t0.1\t0.3")
        self.assertEqual(converted_otu_table[5], "3\t0.0\t0.1\t0.0\t0.0\t-0.4\t-0.4")

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
