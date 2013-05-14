#!/usr/bin/env python

"""Tests of code for performing significance tests of OTU/category associations
"""

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Catherine Lozupone", "Dan Knights", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main
from qiime.otu_category_significance import filter_OTUs, \
    make_contingency_matrix, run_single_G_test, run_G_test_OTUs, \
    add_fdr_correction_to_results, output_results_G_test, \
    run_single_ANOVA, run_ANOVA_OTUs, output_results_ANOVA,\
    run_correlation_OTUs, run_single_correlation, output_results_correlation,\
    get_taxonomy_info, get_category_info,\
    aggregate_multiple_results_ANOVA, aggregate_multiple_results_G_test,\
    aggregate_multiple_results_correlation, get_common_OTUs,\
    test_wrapper_multiple, test_wrapper, get_single_correlation_values,\
    run_paired_T_test_OTUs, run_single_paired_T_test, \
    get_single_paired_T_values, output_results_paired_T_test, sort_rows,\
    sync_mapping_to_otu_table
from numpy import array
from cogent.util.dict2d import Dict2D
from qiime.util import get_tmp_filename
from os import remove
from qiime.parse import parse_mapping_file, parse_otu_table
from biom.parse import parse_biom_table_str
from qiime.format import format_biom_table
from qiime.longitudinal_otu_category_significance import \
        longitudinal_otu_table_conversion_wrapper

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_filter_OTUs(self):
        """filter_OTUs works"""
        otu_table_str = """{"rows": [{"id": "0", "metadata": {}}, {"id": "1",
        "metadata": {}}, {"id": "2", "metadata": {}}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [1, 0, 1.0], [2, 0, 1.0], [2, 1, 1.0], [2, 2, 1.0]], "columns":
        [{"id": "sample1", "metadata": null},
        {"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2541",
        "matrix_type": "sparse", "shape": [3, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T10:13:37.373053", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        #test with all_samples set to True (default), should filter out otu2
        taxonomy_info = get_taxonomy_info(otu_table)
        result = filter_OTUs(otu_table, 0.9)
        self.assertEqual(result, [])
        result = filter_OTUs(otu_table, 0.3333)
        self.assertEqual(result, ['0', '1'])
        
        #test with all_samples set to False. should not filter out otu2
        result = filter_OTUs(otu_table, 0.9, False)
        self.assertEqual(result, ['2'])
        #test with filter set to 0.333 - so should be in 33.3% of samples (1)
        result = filter_OTUs(otu_table, 0.333, False)
        self.assertEqual(result, ['0', '1', '2'])
        #test that is works if a category mapping file is supplied
        cat_mapping = {'sample2': '0', 'sample3': '1'}
        result = filter_OTUs(otu_table, 0.333, category_mapping_info=cat_mapping)
        self.assertEqual(result, ['0'])
        #test that works with a max filter
        result = filter_OTUs(otu_table, 0.333, False, max_filter=0.666667)
        self.assertEqual(result, ['0', '1'])
            
    def test_sync_mapping_to_otu_table(self):
        """sync_mapping_to_otu_table returns a filtered cat mapping"""
        mapping_f4 = """#SampleID\tSomething\tdays_since_epoch
S3\thello\t23
S4\thello\t24
S5\thello\t25
NotInOtuTable1\thello\t26
NotInOtuTable2\thello\t27""".split('\n')
        mapping = parse_mapping_file(mapping_f4)
        
        otu_table_fake2 = """{"rows": [{"id": "0", "metadata": {"Consensus Lineage": "Root;Bacteria"}}, {"id": "3", "metadata": {"Consensus Lineage": "Root;Bacteria;Acidobacteria"}}, {"id": "4", "metadata": {"Consensus Lineage": "Root;Bacteria;Bacteroidetes"}}], "format": "Biological Observation Matrix 1.0.0", "data": [[0, 0, 1.0], [0, 2, 1.0], [1, 0, 2.0], [1, 2, 1.0], [2, 0, 1.0], [2, 2, 9.0]], "columns": [{"id": "S3", "metadata": null}, {"id": "S4", "metadata": null}, {"id": "S5", "metadata": null}], "generated_by": "BIOM-Format 1.0.0-dev", "matrix_type": "sparse", "shape": [3, 3], "format_url": "http://biom-format.org", "date": "2012-07-30T14:36:33.318159", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

        otu_table = parse_biom_table_str(otu_table_fake2)
        new_mapping, removed_samples = sync_mapping_to_otu_table(otu_table,\
                        mapping)
        self.assertEqual(new_mapping[0], [['S3', 'hello', '23'], ['S4', 'hello', '24'], ['S5', 'hello', '25']])
        self.assertEqual(removed_samples, ['NotInOtuTable1', 'NotInOtuTable2']) 

    def test_make_contingency_matrix(self):
        """make_contingency_matrix works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id":
        "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 2, 1.0], [0, 3, 1.0], [1, 0, 1.0], [1, 1, 1.0],
        [2, 2, 1.0], [3, 3, 1.0]], "columns": [{"id": "sample1", "metadata":
        null}, {"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:29:40.964867", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B', 'C']
        result = make_contingency_matrix('0', category_info, otu_table,
                category_values)
        self.assertEqual(result['OTU_pos']['B_pos'], 1)
        self.assertEqual(result['OTU_pos']['C_pos'], 1)
        self.assertEqual(result['OTU_pos']['A_pos'], 0)
        self.assertEqual(result['OTU_neg']['A_pos'], 2)
        self.assertEqual(result['OTU_neg']['B_pos'], 0)
        self.assertEqual(result['OTU_neg']['C_pos'], 0)

    def test_run_single_G_test(self):
        """run_single_G_test works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id":
        "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 2, 1.0], [0, 3, 1.0], [1, 0, 1.0], [1, 1, 1.0],
        [2, 2, 1.0], [3, 3, 1.0]], "columns": [{"id": "sample1", "metadata":
        null}, {"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:29:40.964867", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B', 'C']
        g_val, prob, contingency_matrix = run_single_G_test('0', category_info,
                otu_table, category_values)
        self.assertFloatEqual(g_val, 4.29304060218)
        self.assertFloatEqual(prob, 0.508041627088)
        self.assertEqual(contingency_matrix, {'OTU_pos': {'B_pos': [1, 0.5], 'C_pos': [1, 0.5], 'A_pos': [0, 1.0]}, 'OTU_neg': {'B_pos': [0, 0.5], 'C_pos': [0, 0.5], 'A_pos': [2, 1.0]}})

        #check that it works if samples in mapping are not in OTU table
        otu_table2_str = """{"rows": [{"id": "0", "metadata": null}, {"id":
        "1", "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 1, 1.0], [0, 2, 1.0], [1, 0, 1.0], [2, 1, 1.0],
        [3, 2, 1.0]], "columns": [{"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}, {"id": "sample4", "metadata":
        null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2570",
        "matrix_type": "sparse", "shape": [4, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:36:52.105159", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table2 = parse_biom_table_str(otu_table2_str)

        g_val, prob, contingency_matrix = run_single_G_test('0', category_info,
                otu_table2, category_values)
        self.assertFloatEqual(contingency_matrix, {'OTU_pos': {'B_pos': [1, 0.66666666666666663], 'C_pos': [1, 0.66666666666666663], 'A_pos': [0, 0.66666666666666663]}, 'OTU_neg': {'B_pos': [0, 0.33333333333333331], 'C_pos': [0, 0.33333333333333331], 'A_pos': [1, 0.33333333333333331]}})
        #check that it works is samples in the OTU table are not in the mapping
        category_info2 = {'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        g_val, prob, contingency_matrix = run_single_G_test('0',
                category_info2, otu_table, category_values)
        self.assertFloatEqual(contingency_matrix, {'OTU_pos': {'B_pos': [1, 0.66666666666666663], 'C_pos': [1, 0.66666666666666663], 'A_pos': [0, 0.66666666666666663]}, 'OTU_neg': {'B_pos': [0, 0.33333333333333331], 'C_pos': [0, 0.33333333333333331], 'A_pos': [1, 0.33333333333333331]}})

    def test_run_single_ANOVA(self):
        """run_single_ANOVA works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 5.0], [0, 1, 10.0], [0, 2, 2.0], [0, 3, 1.0],
        [1, 0, 0.0015], [1, 1, 0.0015], [2, 0, 2.0], [2, 2, 1.0], [3, 3, 1.0]],
        "columns": [{"id": "sample1", "metadata": null}, {"id": "sample2",
        "metadata": null}, {"id": "sample3", "metadata": null},
        {"id": "sample4", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2570", "matrix_type": "sparse",
        "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:47:04.545731", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B']
        group_means, prob = run_single_ANOVA('0', category_info, otu_table,
                category_values)
        self.assertEqual(group_means, [7.5, 1.5])
        self.assertFloatEqual(prob, 0.142857142857)
        #test that it works when there are samples in the cat mapping
        #that are not in the OTU table

        otu_table2_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 10.0],
        [0, 1, 2.0], [0, 2, 1.0], [1, 0, 1.0], [2, 1, 1.0], [3, 2, 1.0]],
        "columns": [{"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:49:34.052074", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table2 = parse_biom_table_str(otu_table2_str)

        group_means, prob = run_single_ANOVA('0', category_info,\
            otu_table2, category_values)
        self.assertEqual(group_means, [10.0, 1.5])

        #test that it works when there are samples that are in the OTU
        #table that are not in the category mapping
        category_info2 = {'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}
        group_means, prob = run_single_ANOVA('0', category_info2,\
            otu_table, category_values)
        self.assertEqual(group_means, [10.0, 1.5])
        
    def test_run_single_correlation(self):
        """run_single_correlation works"""
        category_info = {'sample1': '1',
                        'sample2': '2',
                        'sample3': '3',
                        'sample4': '4'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 2.0], [0, 3, 4.0],
        [1, 0, 1.0], [1, 1, 2.0], [1, 2, 3.0], [1, 3, 4.0], [2, 0, 4.0],
        [2, 1, 3.0], [2, 2, 2.0], [2, 3, 1.0], [3, 0, 1.0], [3, 1, 2.0],
        [3, 2, 3.0], [3, 3, 5.0]], "columns": [{"id": "sample1", "metadata":
        null}, {"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:52:44.318249", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B']
        otu_ab_vals, cat_vals = get_single_correlation_values('0', category_info,\
                otu_table)
        r, prob = run_single_correlation(otu_ab_vals, cat_vals)
        self.assertFloatEqual(r, 0.923380516877)
        self.assertFloatEqual(prob, 0.0766194831234)
    
    def test_get_single_correlation_values(self):
        """get_single_correlation_values works
        """
        category_info = {'sample1': '1',
                        'sample2': '2',
                        'sample3': '3',
                        'sample4': '4'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 1.0], [0, 1, 2.0], [0, 2, 2.0], [0, 3, 4.0],
        [1, 0, 1.0], [1, 1, 2.0], [1, 2, 3.0], [1, 3, 4.0], [2, 0, 4.0],
        [2, 1, 3.0], [2, 2, 2.0], [2, 3, 1.0], [3, 0, 1.0], [3, 1, 2.0],
        [3, 2, 3.0], [3, 3, 5.0]], "columns": [{"id": "sample1", "metadata":
        null}, {"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:52:44.318249", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B']
        OTU_abundance_values, category_values = get_single_correlation_values(\
            '0', category_info, otu_table)
        self.assertEqual(OTU_abundance_values, [4.0, 1.0, 2.0, 2.0])
        self.assertEqual(category_values, [4.0, 1.0, 3.0, 2.0])
        #works if the otu table is missing samples

        otu_table2_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null}, {"id":
        "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 2.0],
        [0, 1, 2.0], [0, 2, 4.0], [1, 0, 2.0], [1, 1, 3.0], [1, 2, 4.0],
        [2, 0, 3.0], [2, 1, 2.0], [2, 2, 1.0], [3, 0, 2.0], [3, 1, 3.0],
        [3, 2, 5.0]], "columns": [{"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}, {"id": "sample4", "metadata":
        null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2570",
        "matrix_type": "sparse", "shape": [4, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:56:21.469203", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table2 = parse_biom_table_str(otu_table2_str)

        OTU_abundance_values, category_values = get_single_correlation_values(\
            '0', category_info, otu_table2)
        self.assertEqual(OTU_abundance_values, [4.0, 2.0, 2.0])
        self.assertEqual(category_values, [4.0, 3.0, 2.0])
        #works if the category mapping file is missing samples
        category_info2 = {'sample2': '2',
                        'sample3': '3',
                        'sample4': '4'}
        OTU_abundance_values, category_values = get_single_correlation_values(\
            '0', category_info2, otu_table)
        self.assertEqual(OTU_abundance_values, [4.0, 2.0, 2.0])
        self.assertEqual(category_values, [4.0, 3.0, 2.0])

    def test_get_single_correlation_values_ignore_val(self):
        """get_single_correlation_values works with ignore val"""
        category_info = {'sample1': '1',
                        'sample2': '2',
                        'sample3': '3',
                        'sample4': '4'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 99999.0], [0, 1, 2.0], [0, 2, 2.0], [0, 3, 4.0],
        [1, 0, 1.0], [1, 1, 2.0], [1, 2, 3.0], [1, 3, 4.0], [2, 0, 4.0],
        [2, 1, 3.0], [2, 2, 2.0], [2, 3, 1.0], [3, 0, 1.0], [3, 1, 2.0],
        [3, 2, 3.0], [3, 3, 5.0]], "columns": [{"id": "sample1", "metadata":
        null}, {"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T19:15:22.356507", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        OTU_abundance_values, category_values = \
            get_single_correlation_values('0', category_info, otu_table, \
            ignore_val=99999.0)
        self.assertEqual(OTU_abundance_values, [4.0, 2.0, 2.0])
        self.assertEqual(category_values, [4.0, 3.0, 2.0])

    def test_get_single_paired_T_values(self):
        """get_single_paired_T_values works"""
        cat_mapping = """#SampleID\ttimepoint_zero\tindividual
s1\t1\tA
s2\t0\tA
s3\t1\tB
s4\t0\tB
s5\t1\tC
s6\t0\tC""".split('\n')

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id":
        "1", "metadata": null}, {"id": "2", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 999999999.0],
        [0, 1, 999999999.0], [0, 3, 0.29999999999999999],
        [0, 5, 0.20000000000000001], [1, 1, -0.20000000000000001],
        [1, 2, 999999999.0], [1, 3, 999999999.0], [1, 4, 999999999.0],
        [1, 5, 999999999.0], [2, 1, 0.20000000000000001],
        [2, 3, -0.69999999999999996], [2, 5, 0.10000000000000001]],
        "columns": [{"id": "s1", "metadata": null},
        {"id": "s2", "metadata": null}, {"id": "s3", "metadata": null},
        {"id": "s4", "metadata": null}, {"id": "s5", "metadata": null},
        {"id": "s6", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [3, 6], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:00:27.397644", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)
        mapping_data, header, comments = parse_mapping_file(cat_mapping)
        taxonomy_info = get_taxonomy_info(otu_table)
        OTU_list = ['0', '1', '2']
        
        before_vals, after_vals = get_single_paired_T_values('0', \
            mapping_data, header, 'individual', 'timepoint_zero', \
            otu_table, 999999999.0)
        self.assertFloatEqual(before_vals, [0.0, 0.0])
        self.assertFloatEqual(after_vals, [0.3, 0.2])
        #test of OTU1
        before_vals, after_vals = get_single_paired_T_values('1',
            mapping_data, header, 'individual', 'timepoint_zero',
            otu_table, 999999999.0)
        self.assertFloatEqual(before_vals, [0.0])
        self.assertFloatEqual(after_vals, [-0.2])
        #works when a sample is missing from the OTU table
        #e.g. if an after timepoint dropped out during rarefaction
        #will also drop the before
        otu_table2_str = """{"rows": [{"id": "0", "metadata": null}, {"id":
        "1", "metadata": null}, {"id": "2", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 999999999.0], [0,
        1, 999999999.0], [0, 3, 0.29999999999999999], [1, 1,
        -0.20000000000000001], [1, 2, 999999999.0], [1, 3, 999999999.0], [1, 4,
        999999999.0], [2, 1, 0.20000000000000001], [2, 3,
        -0.69999999999999996]], "columns": [{"id": "s1", "metadata": null},
        {"id": "s2", "metadata": null}, {"id": "s3", "metadata": null}, {"id":
        "s4", "metadata": null}, {"id": "s5", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2564", "matrix_type":
        "sparse", "shape": [3, 5], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:13:14.199954", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""

        otu_table2 = parse_biom_table_str(otu_table2_str)
        before_vals, after_vals = get_single_paired_T_values('0',
            mapping_data, header, 'individual', 'timepoint_zero', otu_table2,
            999999999.0)

        self.assertEqual(before_vals, [0.0])
        self.assertFloatEqual(after_vals, [0.3])
        #works when the before is missing
        otu_table3_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null}],
        "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 999999999.0], [0, 1, 999999999.0],
        [0, 3, 0.29999999999999999], [0, 4, 0.20000000000000001],
        [1, 1, -0.20000000000000001], [1, 2, 999999999.0], [1, 3, 999999999.0],
        [1, 4, 999999999.0], [2, 1, 0.20000000000000001],
        [2, 3, -0.69999999999999996], [2, 4, 0.10000000000000001]], "columns":
        [{"id": "s1", "metadata": null}, {"id": "s2", "metadata": null},
        {"id": "s3", "metadata": null}, {"id": "s4", "metadata": null},
        {"id": "s6", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [3, 5], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:14:52.851951", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table3 = parse_biom_table_str(otu_table3_str)
        before_vals, after_vals = get_single_paired_T_values('0',
            mapping_data, header, 'individual', 'timepoint_zero', otu_table3,
            999999999.0)
        self.assertEqual(before_vals, [0.0])
        self.assertFloatEqual(after_vals, [0.3])
    
    def test_run_single_paired_T_test(self):
        """run_single_paired_T_test works
        """
        cat_mapping = """#SampleID\ttimepoint_zero\tindividual
s1\t1\tA
s2\t0\tA
s3\t1\tB
s4\t0\tB
s5\t1\tC
s6\t0\tC""".split('\n')

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 999999999.0],
        [0, 1, 999999999.0], [0, 3, 0.29999999999999999],
        [0, 5, 0.20000000000000001], [1, 1, -0.20000000000000001],
        [1, 2, 999999999.0], [1, 3, 999999999.0], [1, 4, 999999999.0],
        [1, 5, 999999999.0], [2, 1, 0.20000000000000001],
        [2, 3, -0.69999999999999996], [2, 5, 0.10000000000000001]], "columns":
        [{"id": "s1", "metadata": null}, {"id": "s2", "metadata": null},
        {"id": "s3", "metadata": null}, {"id": "s4", "metadata": null}, {"id":
        "s5", "metadata": null}, {"id": "s6", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2570", "matrix_type":
        "sparse", "shape": [3, 6], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:17:58.541124", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        mapping_data, header, comments = parse_mapping_file(cat_mapping)
        taxonomy_info = get_taxonomy_info(otu_table)
        OTU_list = ['0', '1', '2']
        #should return the results since there should be 66.6% of people to 
        #evaluate  (2 of 3)
        result = run_single_paired_T_test('0', mapping_data, header,
           'individual', 'timepoint_zero', otu_table, 999999999.0, 0.66)
        self.assertEqual(len(result), 4)
        self.assertFloatEqual(result[1], 0.12566591637800242)
        self.assertFloatEqual(result[2], [0.29999999999999999, 0.20000000000000001])
        self.assertEqual(result[3], 2)
        #check the the filter works
        result = run_single_paired_T_test('0', mapping_data, header,
            'individual', 'timepoint_zero', otu_table, 999999999.0, 0.9)
        self.assertEqual(result, None)

    def test_run_paired_T_test_OTUs(self):
        """run_single_paired_T_test works
        """
        cat_mapping = """#SampleID\ttimepoint_zero\tindividual
s1\t1\tA
s2\t0\tA
s3\t1\tB
s4\t0\tB
s5\t1\tC
s6\t0\tC""".split('\n')
        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 999999999.0],
        [0, 1, 999999999.0], [0, 3, 0.29999999999999999],
        [0, 5, 0.20000000000000001], [1, 1, -0.20000000000000001],
        [1, 2, 999999999.0], [1, 3, 999999999.0], [1, 4, 999999999.0],
        [1, 5, 999999999.0], [2, 1, 0.20000000000000001],
        [2, 3, -0.69999999999999996], [2, 5, 0.10000000000000001]], "columns":
        [{"id": "s1", "metadata": null}, {"id": "s2", "metadata": null},
        {"id": "s3", "metadata": null}, {"id": "s4", "metadata": null},
        {"id": "s5", "metadata": null}, {"id": "s6", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2564", "matrix_type":
        "sparse", "shape": [3, 6], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:10:53.571949", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)
        
        mapping_data, header, comments = parse_mapping_file(cat_mapping)
        taxonomy_info = get_taxonomy_info(otu_table)
        OTU_list = ['0', '1', '2']
        all_results = run_paired_T_test_OTUs(OTU_list, mapping_data, header, \
            'individual', 'timepoint_zero', otu_table, 999999999.0, 4)
        self.assertFloatEqual(all_results, {'0': [-4.9999999999999982, 0.12566591637800242, 0.25, 2], '2': [0.46816458878452216, 0.68573031947264562, -0.1333333333333333, 3]})

    def test_output_results_paired_T_test(self):
        """output_results_paired_T_test works
        """
        cat_mapping = """#SampleID\ttimepoint_zero\tindividual
s1\t1\tA
s2\t0\tA
s3\t1\tB
s4\t0\tB
s5\t1\tC
s6\t0\tC""".split('\n')
        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 0, 999999999.0],
        [0, 1, 999999999.0], [0, 3, 0.29999999999999999],
        [0, 5, 0.20000000000000001], [1, 1, -0.20000000000000001],
        [1, 2, 999999999.0], [1, 3, 999999999.0], [1, 4, 999999999.0],
        [1, 5, 999999999.0], [2, 1, 0.20000000000000001],
        [2, 3, -0.69999999999999996], [2, 5, 0.10000000000000001]], "columns":
        [{"id": "s1", "metadata": null}, {"id": "s2", "metadata": null},
        {"id": "s3", "metadata": null}, {"id": "s4", "metadata": null},
        {"id": "s5", "metadata": null}, {"id": "s6", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2564", "matrix_type":
        "sparse", "shape": [3, 6], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:06:21.731328", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)
        
        mapping_data, header, comments = parse_mapping_file(cat_mapping)
        taxonomy_info = get_taxonomy_info(otu_table)
        OTU_list = ['0', '1', '2']
        all_results = run_paired_T_test_OTUs(OTU_list, mapping_data, header, \
            'individual', 'timepoint_zero', otu_table, 999999999.0, 0.6667)
        output = output_results_paired_T_test(all_results)
        self.assertEqual(output, ['OTU\tprob\tT stat\taverage_diff\tnum_pairs\tBonferroni_corrected\tFDR_corrected', '0\t0.125665916378\t-5.0\t0.25\t2\t0.251331832756\t0.251331832756', '2\t0.685730319473\t0.468164588785\t-0.133333333333\t3\t1.37146063895\t0.685730319473'])

    def test_run_ANOVA_OTUs(self):
        """run_ANOVA_OTUs works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}, {"id": "3",
        "metadata": null}], "format": "Biological Observation Matrix v0.9",
        "data": [[0, 0, 5.0], [0, 1, 10.0], [0, 2, 2.0], [0, 3, 1.0],
        [1, 0, 1.0], [1, 3, 2.0], [2, 0, 2.0], [2, 1, 1.0], [2, 2, 10.0],
        [2, 3, 15.0], [3, 0, 1.0], [3, 1, 1.5], [3, 2, 1.3999999999999999],
        [3, 3, 1.3]], "columns": [{"id": "sample1", "metadata": null},
        {"id": "sample2", "metadata": null}, {"id": "sample3", "metadata":
        null}, {"id": "sample4", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:03:44.355433", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B']
        result = run_ANOVA_OTUs(['0', '1', '3'], category_info,\
            otu_table, category_values)
        self.assertEqual(result, {'1': [[0.5, 1.0], 0.69848865542223582, 2.0954659662667074], '0': [[7.5, 1.5], 0.14285714285714285, 0.42857142857142855], '3': [[1.25, 1.3500000000000001], 0.73273875808757438, 2.1982162742627231]})

    def test_output_results_ANOVA(self):
        """output_results_ANOVA works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'B'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}, {"id": "4", "metadata": null}],
        "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 5.0],
        [0, 1, 10.0], [0, 2, 2.0], [0, 3, 1.0], [1, 0, 1.0], [1, 3, 2.0],
        [2, 0, 2.0], [2, 1, 1.0], [2, 2, 10.0], [2, 3, 15.0], [3, 0, 1.0],
        [3, 1, 1.5], [3, 2, 1.3999999999999999], [3, 3, 1.3], [4, 0, 20.0],
        [4, 1, 16.0], [4, 2, 1.3999999999999999], [4, 3, 1.3]], "columns":
        [{"id": "sample1", "metadata": null}, {"id": "sample2", "metadata":
        null}, {"id": "sample3", "metadata": null}, {"id": "sample4",
        "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [5, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T18:00:27.098306", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B']
        taxonomy_info = {'0': 'taxon1',
                        '1': 'taxon2',
                        '2': 'taxon3',
                        '3': 'taxon4',
                        '4': 'taxon5'}
        ANOVA_results = run_ANOVA_OTUs(['0', '2', '1', '3', '4'], category_info,\
            otu_table, category_values)
        output = output_results_ANOVA(ANOVA_results, category_values,\
                        taxonomy_info)
        self.assertEqual(output, ['OTU\tprob\tBonferroni_corrected\tFDR_corrected\tA_mean\tB_mean\tConsensus Lineage', '4\t0.0141325222337\t0.0706626111683\t0.0706626111683\t18.0\t1.35\ttaxon5', '2\t0.0497447318605\t0.248723659303\t0.124361829651\t1.5\t12.5\ttaxon3', '0\t0.142857142857\t0.714285714286\t0.238095238095\t7.5\t1.5\ttaxon1', '1\t0.698488655422\t3.49244327711\t0.873110819278\t0.5\t1.0\ttaxon2', '3\t0.732738758088\t3.66369379044\t0.732738758088\t1.25\t1.35\ttaxon4'])
    
    def test_run_G_test_OTUs(self):
        """run_G_test_OTUs works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 2, 1.0],
        [0, 3, 1.0], [1, 0, 1.0], [1, 1, 1.0], [2, 2, 1.0], [3, 3, 1.0]],
        "columns": [{"id": "sample1", "metadata": null}, {"id": "sample2",
        "metadata": null}, {"id": "sample3", "metadata": null}, {"id":
        "sample4", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:57:12.451685", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B', 'C']
        result = run_G_test_OTUs(['0', '1', '3'], category_info, \
            otu_table, category_values)
        self.assertEqual(result.keys(), ['1', '0', '3'])
        self.assertFloatEqual(result['0'][3], 1.5241248812628101)
        self.assertFloatEqual(result['0'][0], 4.29304060218)
        self.assertFloatEqual(result['0'][1], 0.508041627088)
        self.assertEqual(result['0'][2], {'OTU_pos': {'B_pos': [1, 0.5], 'C_pos': [1, 0.5], 'A_pos': [0, 1.0]}, 'OTU_neg': {'B_pos': [0, 0.5], 'C_pos': [0, 0.5], 'A_pos': [2, 1.0]}})

    def test_fdr_correction_G_test(self):
        """fdr_correction_G_test works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 2, 1.0],
        [0, 3, 1.0], [1, 0, 1.0], [1, 1, 1.0], [1, 2, 1.0], [2, 2, 1.0],
        [2, 3, 1.0], [3, 1, 1.0], [3, 3, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}, {"id": "sample4", "metadata":
        null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2564",
        "matrix_type": "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:52:54.285800", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B', 'C']
        G_test_results = run_G_test_OTUs(['0', '1', '3'], category_info, \
            otu_table, category_values)
        G_test_results = add_fdr_correction_to_results(G_test_results)
        self.assertFloatEqual(G_test_results['0'][-1], 1.52412488126)
        self.assertFloatEqual(G_test_results['1'][-1], 0.938976340277)
        self.assertFloatEqual(G_test_results['3'][-1], 0.828522198394)

    def test_output_results_G_test(self):
        """output_results works"""
        category_info = {'sample1': 'A',
                        'sample2': 'A',
                        'sample3': 'B',
                        'sample4': 'C'}
        otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 2, 1.0],
        [0, 3, 1.0], [1, 0, 1.0], [1, 1, 1.0], [1, 2, 1.0], [2, 2, 1.0],
        [2, 3, 1.0], [3, 1, 1.0], [3, 3, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}, {"id": "sample4", "metadata":
        null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2564",
        "matrix_type": "sparse", "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:52:54.285800", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        category_values = ['A', 'B', 'C']
        taxonomy_info = {'0': 'taxon1',
                        '1': 'taxon2',
                        '2': 'taxon3',
                        '3': 'taxon4'}
        G_test_results = run_G_test_OTUs(['0', '1', '3'], category_info, \
            otu_table, category_values)
        output = output_results_G_test(G_test_results, {})
        self.assertEqual(output, ['OTU\tg_val\tg_prob\tBonferroni_corrected\tFDR_corrected\tOTU_pos##B_pos\tOTU_pos##C_pos\tOTU_pos##A_pos\tOTU_neg##B_pos\tOTU_neg##C_pos\tOTU_neg##A_pos', '0\t4.29304060218\t0.508041627088\t1.52412488126\t1.52412488126\t[1, 0.5]\t[1, 0.5]\t[0, 1.0]\t[0, 0.5]\t[0, 0.5]\t[2, 1.0]', '1\t3.48284992796\t0.625984226851\t1.87795268055\t0.938976340277\t[1, 0.75]\t[0, 0.75]\t[2, 1.5]\t[0, 0.25]\t[1, 0.25]\t[0, 0.5]', '3\t2.14652030109\t0.828522198394\t2.48556659518\t0.828522198394\t[0, 0.5]\t[1, 0.5]\t[1, 1.0]\t[1, 0.5]\t[0, 0.5]\t[1, 1.0]'])
        output = output_results_G_test(G_test_results, taxonomy_info)
        self.assertEqual(output, ['OTU\tg_val\tg_prob\tBonferroni_corrected\tFDR_corrected\tOTU_pos##B_pos\tOTU_pos##C_pos\tOTU_pos##A_pos\tOTU_neg##B_pos\tOTU_neg##C_pos\tOTU_neg##A_pos\tConsensus Lineage', '0\t4.29304060218\t0.508041627088\t1.52412488126\t1.52412488126\t[1, 0.5]\t[1, 0.5]\t[0, 1.0]\t[0, 0.5]\t[0, 0.5]\t[2, 1.0]\ttaxon1', '1\t3.48284992796\t0.625984226851\t1.87795268055\t0.938976340277\t[1, 0.75]\t[0, 0.75]\t[2, 1.5]\t[0, 0.25]\t[1, 0.25]\t[0, 0.5]\ttaxon2', '3\t2.14652030109\t0.828522198394\t2.48556659518\t0.828522198394\t[0, 0.5]\t[1, 0.5]\t[1, 1.0]\t[1, 0.5]\t[0, 0.5]\t[1, 1.0]\ttaxon4'])

    def test_run_correlation_OTUs(self):
        """run_correlation_OTUs works"""
        category_info = {'sample1':'0.1',
                        'sample2':'0.2',
                        'sample3':'0.3',
                        'sample4':'0.4'}

        otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 1.0],
        [0, 2, 2.0], [0, 3, 3.0], [1, 0, 7.0], [1, 1, 5.0], [1, 2, 3.0],
        [1, 3, 1.0], [2, 0, 4.0], [2, 1, 4.2000000000000002], [2, 2, 4.0],
        [2, 3, 4.0], [3, 1, 1.0], [3, 3, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}, {"id": "sample4",
        "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:44:54.905035", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)

        OTU_list = ['1', '0', '2', '3'] 
        result = run_correlation_OTUs(OTU_list, category_info, otu_table)
        #OTU 0 should be positively correlated, 1 negative, and 2&3 neutral
        self.assertFloatEqual(result['0'][0], 1.0)#r
        self.assertFloatEqual(result['0'][1], 0.0)#prob
        self.assertEqual(result['0'][2], "[3.0, 0.0, 2.0, 1.0]")
        x = result['0'][3]
        x = x[1:-1]
        x = x.split(',')
        x = [float(i) for i in x]
        self.assertFloatEqual(x, [0.40000000000000002, 0.10000000000000001, 0.29999999999999999, 0.20000000000000001])
        self.assertFloatEqual(result['1'][0], -0.99999999999999956)
        self.assertFloatEqual(result['1'][1], 4.4408920985006281e-16)
        y = result['1'][2]
        y = y[1:-1]
        y = y.split(',')
        y = [float(i) for i in y]
        self.assertFloatEqual(y, [1.0, 7.0, 3.0, 5.0])

        #test that appropriate error is raised is categorical
        category_info = {'sample1':'A',
                        'sample2':'B',
                        'sample3':'A',
                        'sample4':'B'}
        self.assertRaises(ValueError, run_correlation_OTUs, OTU_list,
                          category_info, otu_table)
        
    def test_output_results_correlation(self):
        """output_results_correlation works"""
        category_info = {'sample1':'0.1',
                        'sample2':'0.2',
                        'sample3':'0.3',
                        'sample4':'0.4'}
        otu_table_str = """{"rows": [{"id": "0", "metadata": null},
        {"id": "1", "metadata": null}, {"id": "2", "metadata": null},
        {"id": "3", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 1.0],
        [0, 2, 2.0], [0, 3, 3.0], [1, 0, 7.0], [1, 1, 5.0], [1, 2, 3.0],
        [1, 3, 1.0], [2, 0, 4.0], [2, 1, 4.2000000000000002], [2, 2, 4.0],
        [2, 3, 4.0], [3, 1, 1.0], [3, 3, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}, {"id": "sample4",
        "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [4, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:44:54.905035", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        otu_table = parse_biom_table_str(otu_table_str)
        
        taxonomy_info = {'0': 'taxon1',
                        '1': 'taxon2',
                        '2': 'taxon3',
                        '3': 'taxon4'}
        OTU_list = ['1', '0', '2', '3'] 
        result = run_correlation_OTUs(OTU_list, category_info, otu_table)
        output = output_results_correlation(result, taxonomy_info)
        self.assertEqual(output[0], 'OTU\tprob\totu_values_y\tcat_values_x\tBonferroni_corrected\tFDR_corrected\tr\tConsensus Lineage')
        line2 = output[2].split('\t')
        self.assertEqual(line2[0], '1')
        self.assertFloatEqual(float(line2[1]), 4.4408920985e-16)
        self.assertFloatEqual(float(line2[4]), 1.7763568394e-15)
        self.assertFloatEqual(float(line2[5]), 8.881784197e-16)
        self.assertFloatEqual(float(line2[6]), -1.0)
        self.assertFloatEqual(line2[7], 'taxon2')
        
        self.assertEqual(len(output), 5)

    def test_get_taxonomy_info(self):
        """get_taxonomy_info works"""
        otu_table_str = """{"rows": [{"id": "0", "metadata": null}, {"id": "1",
        "metadata": null}, {"id": "2", "metadata": null}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [1, 0, 1.0], [2, 0, 1.0], [2, 1, 1.0], [2, 2, 1.0]], "columns":
        [{"id": "sample1", "metadata": null}, {"id": "sample2", "metadata":
        null}, {"id": "sample3", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [3, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:27:12.329308", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        
        otu_table = parse_biom_table_str(otu_table_str)
        taxonomy_info = get_taxonomy_info(otu_table)
        self.assertEqual(taxonomy_info, {})

        #test that it parses otu tables with taxonomy fields appropriately
        otu_table_str = """{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["Bacteria", "Bacteroidetes", "Bacteroidales", "Parabacteroidaceae",
        "Unclassified", "otu_475"]}}, {"id": "1", "metadata": {"taxonomy":
        ["Bacteria", "Bacteroidetes", "Bacteroidales", "adhufec77-25",
        "Barnesiella", "Barnesiella_viscericola", "otu_369"]}}, {"id": "2",
        "metadata": {"taxonomy": ["Bacteria", "Firmicutes", "Clostridia",
        "Clostridiales", "Faecalibacterium", "Unclassified", "otu_1121"]}}],
        "format": "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [1, 0, 1.0], [2, 0, 1.0], [2, 1, 1.0], [2, 2, 1.0]], "columns":
        [{"id": "sample1", "metadata": null}, {"id": "sample2", "metadata":
        null}, {"id": "sample3", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [3, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:31:44.658057", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}"""
        
        otu_table = parse_biom_table_str(otu_table_str)
        taxonomy_info = get_taxonomy_info(otu_table)
        self.assertEqual(taxonomy_info, {'1':
            'Bacteria; Bacteroidetes; Bacteroidales; adhufec77-25; ' + \
            'Barnesiella; Barnesiella_viscericola; otu_369', '0':
            'Bacteria; Bacteroidetes; Bacteroidales; Parabacteroidaceae; ' + \
            'Unclassified; otu_475', '2':
            'Bacteria; Firmicutes; Clostridia; Clostridiales; ' + \
            'Faecalibacterium; Unclassified; otu_1121'})

    def test_get_category_info(self):
        """get_category_info works"""
        category_mapping = """#SampleID\tcat1\tcat2
sample1\tA\t0
sample2\tB\t8.0
sample3\tC\t1.0""".split('\n')
        mapping_data, header, comments = parse_mapping_file(category_mapping)
        result, cat_vals = get_category_info(mapping_data, header, 'cat1')
        self.assertEqual(result, {'sample1': 'A', 'sample3': 'C', 'sample2': 'B'})
        self.assertEqual(cat_vals, (['A', 'B', 'C']))
        mapping_data, header, comments = parse_mapping_file(category_mapping)
        result, cat_vals = get_category_info(mapping_data, header, \
                        'cat2', threshold=5.0)
        self.assertEqual(result, {'sample1': '0', 'sample3': '0', 'sample2': '1'})
        self.assertEqual(cat_vals, (['0', '1']))
        
    def test_test_wrapper(self):
        """runs the specified statistical test"""
        otu_table1 = parse_biom_table_str("""{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["lineage0"]}}, {"id": "1", "metadata": {"taxonomy": ["lineage1"]}},
        {"id": "2", "metadata": {"taxonomy": ["lineage2"]}}],
        "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0],
        [0, 1, 2.0], [0, 3, 1.0], [1, 0, 1.0], [1, 3, 1.0], [2, 0, 1.0],
        [2, 1, 1.0], [2, 2, 1.0], [2, 3, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null}, {"id":
        "sample3", "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2556", "matrix_type":
        "sparse", "shape": [3, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T15:42:03.286885", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}""")
        
        category_mapping = ['#SampleID\tcat1\tcat2',
                                      'sample1\tA\t0',
                                      'sample2\tA\t8.0',
                                      'sample3\tB\t1.0',
                                      'sample4\tB\t1.0']
        category_mapping = parse_mapping_file(category_mapping)
        OTU_list = ['1', '0'] 

        # ANOVA
        threshold = None
        filt = 0
        otu_include = None
        category = 'cat1'

        # get expected ANOVA output from file
        results1 = test_wrapper('ANOVA', otu_table1, category_mapping, \
            category, threshold, filt, otu_include=None)
        self.assertEqual(len(results1), 4)
        self.assertEqual(results1[0], 'OTU\tprob\tBonferroni_corrected\tFDR_corrected\tA_mean\tB_mean\tConsensus Lineage')
        otu0_results = results1[1].split('\t')
        A_mean = float(otu0_results[4])
        self.assertFloatEqual(A_mean, 0.5)
        B_mean = float(otu0_results[5])
        self.assertFloatEqual(B_mean, 0.166666666666666)


        # get expected ANOVA means when it is specified not to convert the table to relative abundance
        results1 = test_wrapper('ANOVA', otu_table1, category_mapping, \
            category, threshold, filt, otu_include=None, otu_table_relative_abundance=True)
        self.assertEqual(len(results1), 4)
        self.assertEqual(results1[0], 'OTU\tprob\tBonferroni_corrected\tFDR_corrected\tA_mean\tB_mean\tConsensus Lineage')
        otu0_results = results1[1].split('\t')
        A_mean = float(otu0_results[4])
        self.assertFloatEqual(A_mean, 1.5)
        B_mean = float(otu0_results[5])
        self.assertFloatEqual(B_mean, 0.5)
        
        # correlation
        threshold = None
        otu_include = None
        results1 = test_wrapper('correlation', otu_table1, category_mapping, \
            'cat2', threshold, filt=0, otu_include=None)
        self.assertEqual(len(results1), 4)
        self.assertEqual(results1[0], 'OTU\tprob\totu_values_y\tcat_values_x\tBonferroni_corrected\tFDR_corrected\tr\tConsensus Lineage')
        
        #correlation with filter: should filter out 1 of the 3 otus when
        #set to >0.75
        filt = 0.75
        results1 = test_wrapper('correlation', otu_table1, category_mapping, \
            'cat2', threshold, filt=filt)
        self.assertEqual(len(results1), 3)
        #should filter out 2 of 3 when sent to 0.9
        results1 = test_wrapper('correlation', otu_table1, category_mapping, \
            'cat2', threshold, filt=0.9)
        self.assertEqual(len(results1), 2)
        #should filter out none when sent to 0.5
        results1 = test_wrapper('correlation', otu_table1, category_mapping, \
            'cat2', threshold, filt=0.5)
        self.assertEqual(len(results1), 4)

    def test_aggregate_multiple_results_ANOVA(self):
        """aggregate_multiple_results_ANOVA works"""
        all_results = {
            '0':[[[10,20],.05,.1], [[12,18],.1,.2]],
            '1':[[[15,14],.5,1], [[12,11],.6, 1.2]]
            }
        result = aggregate_multiple_results_ANOVA(all_results)
        self.assertFloatEqual(result['0'], [array([11,19]),.075, .15])
        self.assertFloatEqual(result['1'], [array([13.5,12.5]),.55, 1.1])
                       
    def test_aggregate_multiple_results_G_test(self):
        """aggregate_multiple_results_G_test works"""
        cm_0_0 = Dict2D({'OTU_pos':{'A_pos':2, 'B_pos':0},
                  'OTU_neg':{'A_pos':3,'B_pos':4}},Default=0, Pad=True)
        cm_0_1 = Dict2D({'OTU_pos':{'A_pos':5, 'B_pos':1},
                  'OTU_neg':{'A_pos':6,'B_pos':10}},Default=0, Pad=True)
        cm_1_0 = Dict2D({'OTU_pos':{'A_pos':3, 'B_pos':5},
                  'OTU_neg':{'A_pos':1,'B_pos':1}},Default=0, Pad=True)
        cm_1_1 = Dict2D({'OTU_pos':{'A_pos':4, 'B_pos':4},
                  'OTU_neg':{'A_pos':0,'B_pos':1}},Default=0, Pad=True)

        exp_cm_0 = Dict2D({'OTU_pos':{'A_pos':3.5, 'B_pos':0.5},
                  'OTU_neg':{'A_pos':4.5,'B_pos':7}},Default=0, Pad=True)
        exp_cm_1 = Dict2D({'OTU_pos':{'A_pos':3.5, 'B_pos':4.5},
                  'OTU_neg':{'A_pos':.5,'B_pos':1}},Default=0, Pad=True)
        all_results = {\
                 '0':[[.5, .05, cm_0_0, .1], [.6,.1,cm_0_1,.2]],
                 '1':[[.6, .5, cm_1_0, 1], [.8,.6, cm_1_1, 1.2]]
                 }
        results = aggregate_multiple_results_G_test(all_results)

        self.assertFloatEqual(results['0'], [0.55, 0.075, exp_cm_0, .15])
        self.assertFloatEqual(results['1'], [0.7, 0.55, exp_cm_1, 1.1])
        

    def test_aggregate_multiple_results_correlation(self):
        """aggregate_multiple_results_correlation works"""
        all_results = {
            '0':[[-1,.05,'NA', 'NA',.1], [-.5,.1,'NA','NA',.2]],
            '1':[[.4,.5,'NA','NA',1], [.5,.6,'NA','NA', 1.2]]
            }
        result = aggregate_multiple_results_correlation(all_results)
        self.assertFloatEqual(result['0'], [-.75, .075, 'NA', 'NA'])
        self.assertFloatEqual(result['1'], [.45,.55, 'NA', 'NA'])


    def test_get_common_OTUs(self):
        """get_common_OTUs works"""

        # create the temporary OTU tables
        otu_table1 = parse_biom_table_str("""{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["lineage0"]}}, {"id": "1", "metadata": {"taxonomy": ["lineage1"]}},
        {"id": "2", "metadata": {"taxonomy": ["lineage2"]}}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [1, 0, 1.0], [2, 0, 1.0], [2, 1, 1.0], [2, 2, 1.0]], "columns":
        [{"id": "sample1", "metadata": null}, {"id": "sample2", "metadata":
        null}, {"id": "sample3", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [3, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:19:22.290776", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}""")
        
        otu_table2 = parse_biom_table_str("""{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["lineage0"]}}, {"id": "1", "metadata": {"taxonomy": ["lineage1"]}},
        {"id": "2", "metadata": {"taxonomy": ["lineage2"]}}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [1, 0, 1.0], [2, 1, 1.0], [2, 2, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null},
        {"id": "sample3", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [3, 3], "format_url":
        "http://www.qiime.org/svn_documentation/biom_format.html",
        "date": "2011-12-21T17:21:13.685763", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}""")
        
        otu_table3 = parse_biom_table_str("""{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["lineage0"]}}, {"id": "2", "metadata": {"taxonomy": ["lineage2"]}}],
        "format": "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [1, 0, 1.0], [1, 1, 1.0], [1, 2, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null}, {"id":
        "sample3", "metadata": null}], "generated_by":
        "QIIME 1.4.0-dev, svn revision 2564", "matrix_type": "sparse",
        "shape": [2, 3], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T17:22:55.602729", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}""")
        
        category_info = {'sample1':'0.1',
                        'sample2':'0.2',
                        'sample3':'0.3'}
        OTU_list = ['1', '0', '2'] 


        # case where one OTU is missing from one file
        otu_tables = [otu_table1, otu_table2, otu_table3]
        filt = 0
        filter_all_samples = False
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_tables, filt, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['0','2'])
        exp_taxonomy = {'0': 'lineage0', '2': 'lineage2'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

        # case where no OTUs should be filtered
        otu_tables = [otu_table1, otu_table2]
        filt = 0
        filter_all_samples = False
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_tables, filt, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['0','1','2'])
        exp_taxonomy = {'0': 'lineage0', '1':'lineage1', '2': 'lineage2'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

        # case where one OTU should be filtered due to filter_all_samples
        otu_tables = [otu_table1,otu_table2]
        filt = 0
        filter_all_samples = True
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_tables, filt, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['0','1'])
        exp_taxonomy = {'0': 'lineage0', '1':'lineage1'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

        # case where two OTUs should be filtered due to filt value
        otu_tables = [otu_table1,otu_table2]
        filt = 0.66666667
        filter_all_samples = False
        otu_include = None
        OTU_list, taxonomy = get_common_OTUs(otu_tables, filt, \
                                             category_info, \
                                             filter_all_samples, \
                                             otu_include=otu_include)
        exp_OTU_list = sorted(['2'])
        exp_taxonomy = {'2':'lineage2'}
        self.assertEqual(sorted(OTU_list), exp_OTU_list)
        self.assertEqual(taxonomy,exp_taxonomy)

    def test_test_wrapper_multiple(self):
        """test_wrapper_multiple works"""
        # create the temporary OTU tables
        otu_table1 = parse_biom_table_str("""{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["lineage0"]}}, {"id": "1", "metadata": {"taxonomy": ["lineage1"]}},
        {"id": "2", "metadata": {"taxonomy": ["lineage2"]}}],
        "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0],
        [0, 1, 2.0], [0, 3, 1.0], [1, 0, 1.0], [1, 3, 1.0], [2, 0, 1.0],
        [2, 1, 1.0], [2, 2, 1.0], [2, 3, 1.0]], "columns": [{"id": "sample1",
        "metadata": null}, {"id": "sample2", "metadata": null}, {"id":
        "sample3", "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2556", "matrix_type":
        "sparse", "shape": [3, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T15:42:03.286885", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}""")
        
        otu_table2 = parse_biom_table_str("""{"rows": [{"id": "0", "metadata": {"taxonomy":
        ["lineage0"]}}, {"id": "1", "metadata": {"taxonomy": ["lineage1"]}},
        {"id": "2", "metadata": {"taxonomy": ["lineage2"]}}], "format":
        "Biological Observation Matrix v0.9", "data": [[0, 1, 2.0],
        [0, 3, 1.0], [1, 0, 1.0], [1, 3, 1.0], [2, 1, 1.0], [2, 2, 1.0],
        [2, 3, 1.0]], "columns": [{"id": "sample1", "metadata": null},
        {"id": "sample2", "metadata": null}, {"id": "sample3",
        "metadata": null}, {"id": "sample4", "metadata": null}],
        "generated_by": "QIIME 1.4.0-dev, svn revision 2556", "matrix_type":
        "sparse", "shape": [3, 4], "format_url":
        "http://www.qiime.org/svn_documentation/documentation/biom_format.html",
        "date": "2011-12-21T15:45:47.278634", "type": "OTU table", "id": null,
        "matrix_element_type": "float"}""")

        category_mapping = ['#SampleID\tcat1\tcat2',
                                      'sample1\tA\t0',
                                      'sample2\tA\t8.0',
                                      'sample3\tB\t1.0',
                                      'sample4\tB\t1.0']
        category_mapping = parse_mapping_file(category_mapping)
        OTU_list = ['1', '0'] 

        # ANOVA
        otu_tables = [otu_table1,otu_table2]
        threshold = None
        filt = 0
        otu_include = None
        category = 'cat1'

        # get expected ANOVA output from each file separately
        results1 = test_wrapper('ANOVA', otu_table1, category_mapping, category, threshold, \
                 filt, otu_include=None)
        results2 = test_wrapper('ANOVA', otu_table2, category_mapping, category, threshold, \
                 filt, otu_include=None)

        results = test_wrapper_multiple('ANOVA', otu_tables,
                                        category_mapping, category,
                                        threshold, filt,
                                        otu_include)

        # multiple results should be combination of individual results
        results1 = [result.split('\t') for result in results1]        
        results1_dict = {}
        for line in results1[1:]:
            results1_dict[line[0]] = line[1:]
        results2 = [result.split('\t') for result in results2]
        results2_dict = {}
        for line in results2[1:]:
            results2_dict[line[0]] = line[1:]
        results = [result.split('\t') for result in results]
        results_dict = {}
        for line in results[1:]:
            results_dict[line[0]] = line[1:]
    
        # for every OTU in the combined results, compare to average of separate
        for k in results_dict.keys():
            # skip FDR corrected values because they are calculated _after_ the
            # method combines separate results
            for i in [0,1,3,4]:
                entry1 = float(results1_dict[k][i])
                entry2 = float(results2_dict[k][i])
                entry_combined = float(results_dict[k][i])
                mean = round((entry1+entry2)/2.0, 3)
                self.assertEqual(round(entry_combined, 3), mean)

        # correlation
        threshold = None
        filt = 0
        otu_include = None
        category = 'cat2'

        # get expected correlation output from each file separately
        results1 = test_wrapper('correlation', otu_table1, category_mapping, category, threshold, \
                 filt, otu_include=None)
        results2 = test_wrapper('correlation', otu_table2, category_mapping, category, threshold, \
                 filt, otu_include=None)

        results = test_wrapper_multiple('correlation', otu_tables,
                                        category_mapping, category,
                                        threshold, filt,
                                        otu_include)

        # multiple results should be combination of individual results
        results1 = [result.split('\t') for result in results1]        
        results1_dict = {}
        for line in results1[1:]:
            results1_dict[line[0]] = line[1:]
        results2 = [result.split('\t') for result in results2]
        results2_dict = {}
        for line in results2[1:]:
            results2_dict[line[0]] = line[1:]
        results = [result.split('\t') for result in results]
        results_dict = {}
        for line in results[1:]:
            results_dict[line[0]] = line[1:]
    
        # for every OTU in the combined results, compare to average of separate
        for k in results_dict.keys():
            # skip FDR corrected values and bonferroni values because they are 
            #calculated _after_ the
            # method combines separate results
            for i in [0,]:
                entry1 = float(results1_dict[k][i])
                entry2 = float(results2_dict[k][i])
                entry_combined = float(results_dict[k][i])
                mean = round((entry1+entry2)/2.0, 3)
                self.assertEqual(round(entry_combined, 3), mean)


        # G test
        threshold = None
        filt = 0
        otu_include = None
        category = 'cat1'

        # get expected G_TEST output from each file separately
        results1 = test_wrapper('g_test', otu_table1, category_mapping, category, threshold, \
                 filt, otu_include=['0'])

        results2 = test_wrapper('g_test', otu_table2, category_mapping, category, threshold, \
                 filt, otu_include=['0'])
        results = test_wrapper_multiple('g_test', otu_tables,
                                        category_mapping, category,
                                        threshold, filt,
                                        otu_include=['0'])

        # multiple results should be combination of individual results
        results1 = [result.split('\t') for result in results1]        
        results1_dict = {}
        for line in results1[1:]:
            results1_dict[line[0]] = line[1:]
        results2 = [result.split('\t') for result in results2]
        results2_dict = {}
        for line in results2[1:]:
            results2_dict[line[0]] = line[1:]
        results = [result.split('\t') for result in results]
        results_dict = {}
        for line in results[1:]:
            results_dict[line[0]] = line[1:]
    
        # convert all numeric entries to floats rounded to 5 decimal places
        for k in results_dict.keys():
            # skip FDR corrected values because they are calculated _after_ the
            # method combines separate results
            for i in [0,1,2]:
                entry1 = float(results1_dict[k][i])
                entry2 = float(results2_dict[k][i])
                entry_combined = float(results_dict[k][i])
                mean = round((entry1+entry2)/2.0, 3)

    def test_sort_rows(self):
        """sort_rows works"""
        output = ['OTU\tprob\tBonferroni_corrected\tFDR_corrected\tA_mean\tB_mean\tConsensus Lineage', '1\t0.698488655422\t3.49244327711\t0.873110819278\t0.5\t1.0\ttaxon2', '0\t0.142857142857\t0.714285714286\t0.238095238095\t7.5\t1.5\ttaxon1', '3\t0.732738758088\t3.66369379044\t0.732738758088\t1.25\t1.35\ttaxon4', '2\t0.0497447318605\t0.248723659303\t0.124361829651\t1.5\t12.5\ttaxon3', '4\t0.0141325222337\t0.0706626111683\t0.0706626111683\t18.0\t1.35\ttaxon5']
        result = sort_rows(output, 1)
        self.assertEqual(result, ['OTU\tprob\tBonferroni_corrected\tFDR_corrected\tA_mean\tB_mean\tConsensus Lineage', '4\t0.0141325222337\t0.0706626111683\t0.0706626111683\t18.0\t1.35\ttaxon5', '2\t0.0497447318605\t0.248723659303\t0.124361829651\t1.5\t12.5\ttaxon3', '0\t0.142857142857\t0.714285714286\t0.238095238095\t7.5\t1.5\ttaxon1', '1\t0.698488655422\t3.49244327711\t0.873110819278\t0.5\t1.0\ttaxon2', '3\t0.732738758088\t3.66369379044\t0.732738758088\t1.25\t1.35\ttaxon4'])

    def test_longitudinal_correlation_filter(self):
        """longitudinal correlation filter works as expected
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
        converted_otu_table = longitudinal_otu_table_conversion_wrapper(\
            otu_table, category_mapping, 'individual', 'timepoint_zero')
        output = test_wrapper('correlation', converted_otu_table,\
            category_mapping, category='timepoint', threshold=None, filt=0.1,\
            ignore_val=999999999.0, otu_table_relative_abundance=True)
        self.assertEqual(len(output), 5)
        output = test_wrapper('correlation', converted_otu_table,\
            category_mapping, category='timepoint', threshold=None, filt=0.6,\
            ignore_val=999999999.0, otu_table_relative_abundance=True)
        self.assertEqual(len(output), 3)


#run unit tests if run from command-line
if __name__ == '__main__':
    main()
