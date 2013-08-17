#!/usr/bin/env python
# File created on 16 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main 
from qiime.ocs import sync_biom_and_mf, get_sample_cats, get_cat_sample_groups, \
    get_sample_indices, row_generator, run_ocs_test, fdr_correction, \
    bonferroni_correction, output_formatter, sort_by_pval
from qiime.pycogent_backports.test import (parametric_correlation_significance,
    nonparametric_correlation_significance, fisher_confidence_intervals,
    pearson, spearman, kendall_correlation, G_fit, ANOVA_one_way, 
    kruskal_wallis, mw_test, mw_boot)
from cogent.maths.stats.test import (t_two_sample, mc_t_two_sample)
from numpy import array
from numpy.random import seed
from cogent.util.dict2d import Dict2D
from qiime.util import get_tmp_filename
from os import remove
from qiime.parse import parse_mapping_file_to_dict, parse_otu_table
from biom.parse import parse_biom_table_str
from qiime.format import format_biom_table

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def test_sync_biom_and_mf(self):
        """sync_biom_and_mf works"""

        #set up otu table and mapping file 
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
        "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
        "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
        [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
        [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
        [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
        [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
        [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
        "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
        "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
        ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
        "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
        {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
        {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
        {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
        {"id": "Sample6", "metadata": null}]}"""

        mapping_file_1 = """#SampleID\ttest_cat\ttest_corr
        Sample1\tcat1\t1
        Sample2\tcat1\t2
        Sample3\tcat2\t3
        Sample4\tcat2\t4
        Sample5\tcat3\t5
        Sample6\tcat3\t6
        NotInOtuTable1\tcat5\t7
        NotInOtuTable2\tcat5\t8""".split('\n')
        
        mf, _ = parse_mapping_file_to_dict(mapping_file_1)
        bt = parse_biom_table_str(otu_table_1)

        pmf_out = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
        'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
        'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
        'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
        'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
        'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}

        bt_out = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "QIIME 1.7.0-dev, ocs_collab@c4c972f","date": "2013-08-16T11:14:06.887109","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}"""
        bt_out = parse_biom_table_str(bt_out)

        npmf_result, nbt_result = sync_biom_and_mf(mf, bt)
        self.assertEqual(npmf_result, pmf_out)
        self.assertEqual(nbt_result, bt_out)

    def test_get_sample_cats(self):
        """get_sample_cats works"""
        pmf_in = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
            'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
            'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
            'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
            'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
            'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}

        get_sample_cats_out = {'Sample1': 'cat1',
        'Sample2': 'cat1',
        'Sample3': 'cat2',
        'Sample4': 'cat2',
        'Sample5': 'cat3',
        'Sample6': 'cat3'}

        get_sample_cats_result = get_sample_cats(pmf_in, 'test_cat')
        self.assertEqual(get_sample_cats_result, get_sample_cats_out)

    def test_get_cat_sample_groups(self):
        """get_cat_sample_groups works"""

        sample_cats = {'Sample1': 'cat1',
        'Sample2': 'cat1',
        'Sample3': 'cat2',
        'Sample4': 'cat2',
        'Sample5': 'cat3',
        'Sample6': 'cat3'}

        cat_sample_groups_out = {'cat1': ['Sample1', 'Sample2'],
            'cat2': ['Sample4', 'Sample3'],
            'cat3': ['Sample5', 'Sample6']}

        cat_sample_groups_result = get_cat_sample_groups(sample_cats)
        self.assertEqual(cat_sample_groups_result, cat_sample_groups_out)

    def test_get_sample_indices(self):
        """get_sample_indices works"""

        cat_sample_groups = {'cat1': ['Sample1', 'Sample2'],
            'cat2': ['Sample4', 'Sample3'],
            'cat3': ['Sample5', 'Sample6']}

        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
        "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
        "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
        [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
        [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
        [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
        [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
        [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
        "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
        "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
        ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
        "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
        {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
        {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
        {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
        {"id": "Sample6", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)

        sample_indices_out = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        sample_indices_result = get_sample_indices(cat_sample_groups, bt)
        self.assertEqual(sample_indices_result, sample_indices_out)

    def test_row_generator(self):
        """row_generator works"""

        sample_indices = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
        "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
        "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
        [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
        [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
        [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
        [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
        [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
        "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
        "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
        ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
        "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
        {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
        {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
        {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
        {"id": "Sample6", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)

        all_rows_result = []
        row_data_gen = row_generator(bt, sample_indices)
        for row in row_data_gen:
            all_rows_result.append(row)

        row_0_out = [array([ 28.,  52.]), array([ 78.,  51.]), array([ 16.,  77.])]
        row_1_out = [array([ 25.,  14.]), array([ 32.,  11.]), array([ 48.,  63.])]
        self.assertEqual(all_rows_result[0], row_0_out)
        self.assertEqual(all_rows_result[1], row_1_out)

    def test_run_ocs_test_nonparametric_t_test(self):
        """run_ocs_test_nonparametric_t_test works"""
        
        test_choices = {'nonparametric_t_test': mc_t_two_sample}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
        "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
        "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
        [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
        [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
        [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
        [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
        [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
        "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
        "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
        ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
        "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
        {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
        {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
        {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
        {"id": "Sample6", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        row_gen = row_generator(bt, sample_indices)

        test_stats_output = [-1.3564095042777529,
            -0.16872981624383912,
            -0.83205029433784361,
            -0.95713665769637957,
            -1.3375524861268122,
            -2.2346432765093422]

        pvals_output = [0.66600000000000004,
         1.0,
         0.68300000000000005,
         0.68000000000000005,
         0.65400000000000003,
         0.32000000000000001]

        means_output = [[40.0, 64.5, 46.5],
             [19.5, 21.5, 55.5],
             [16.5, 42.0, 45.5],
             [52.0, 67.5, 47.5],
             [28.5, 49.5, 9.0],
             [20.0, 76.0, 39.5]]

        seed(0) # set rand generator so that results are comparable
        test_stats_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'nonparametric_t_test', test_choices)

        self.assertFloatEqual(test_stats_result, test_stats_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    def test_run_ocs_test_boodstrap_mann_whitney_u(self):
        """run_ocs_test_boodstrap_mann_whitney_u works"""

        test_choices = {'bootstrap_mann_whitney_u': mw_boot}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = row_generator(bt, sample_indices)   

        test_stat_output = [10.0, 15.0, 11.0, 14.0, 15.0, 9.0]
        pvals_output = [0.603, 0.029, 0.4, 0.089, 0.034, 0.82]
        means_output = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]

        seed(0) # set rand generator so that results are comparable
        test_stat_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'bootstrap_mann_whitney_u', test_choices)

        self.assertFloatEqual(test_stat_result, test_stat_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    def test_run_ocs_test_parametric_t_test(self):
        """run_ocs_test_parametric_t_test works"""
        
        test_choices = {'parametric_t_test': t_two_sample}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = row_generator(bt, sample_indices)

        test_stat_output = [0.43577690622483684,
             -2.5911938781738648,
             -1.3573515147239095,
             1.2101173913086851,
             2.137178815882979,
             0.0099191576638653078]
        pvals_output = [0.67823972846362579,
             0.041145883121579255,
             0.2235024418313547,
             0.27174025956151748,
             0.076447615888438444,
             0.9924073718332862]
        means_output = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]

        test_stat_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'parametric_t_test', test_choices)

        self.assertFloatEqual(test_stat_result, test_stat_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    def test_run_ocs_test_mann_whitney_u(self):
        """run_ocs_test_boodstrap_mann_whitney_u works"""

        test_choices = {'mann_whitney_u': mw_test}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = row_generator(bt, sample_indices)

        test_stat_output = [10.0, 15.0, 11.0, 14.0, 15.0, 9.0]
        pvals_output = [0.5637028616507731,
             0.043308142810791955,
             0.38363032713198975,
             0.083264516663550406,
             0.043308142810791955,
             0.77282999268444752]
        means_output = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]
        
        test_stat_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'mann_whitney_u', test_choices)
        
        self.assertFloatEqual(test_stat_result, test_stat_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    # ANOVA, G_Fit, kruskal_wallis
    def test_run_ocs_test_ANOVA(self):
        """run_ocs_test_ANOVA works"""
        test_choices = {'ANOVA': ANOVA_one_way}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = row_generator(bt, sample_indices)

        test_stat_output = [0.18990151199889027,
             6.7142857142857144,
             1.8424031345232912,
             1.4643841007477372,
             4.5675332910589734,
             9.8389688760617899e-05]
        pvals_output = [0.6782397284636259,
             0.041145883121579234,
             0.22350244183135481,
             0.27174025956151771,
             0.076447615888438403,
             0.9924073718332751]
        means_output = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]
        
        test_stat_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'ANOVA', test_choices)
        
        self.assertFloatEqual(test_stat_result, test_stat_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    def test_run_ocs_test_g_fit(self):
        """run_ocs_test_g_fit works"""
        test_choices = {'g_test': G_fit}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = row_generator(bt, sample_indices)

        test_stat_output = [0.8950130401309585,
             8.6948783805472942,
             6.5397009199496443,
             2.2281537448054953,
             11.541070115516771,
             0.00064935138712822981]
        pvals_output = [0.34412242732851783,
             0.0031910540870178925,
             0.010549308294222293,
             0.13551569348660794,
             0.00068075444949030543,
             0.97967020739471489]
        means_output = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]

        test_stat_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'g_test', test_choices)
        
        self.assertFloatEqual(test_stat_result, test_stat_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    def test_run_ocs_test_kruskal_wallis(self):
        """run_ocs_test_kruskal_wallis works"""
        test_choices = {'kruskal_wallis': kruskal_wallis}
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

        bt = parse_biom_table_str(otu_table_1)
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = row_generator(bt, sample_indices)

        test_stat_output = [0.33333333333333215,
             4.0833333333333321,
             0.75903614457831325,
             3.0,
             4.0833333333333321,
             0.083333333333332149]
        pvals_output = [0.56370286165077377,
             0.043308142810792101,
             0.38363032713198986,
             0.08326451666355042,
             0.043308142810792101,
             0.77282999268444919]
        means_output = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]
        
        test_stat_result, pvals_result, means_result = \
            run_ocs_test(row_gen, 'kruskal_wallis', test_choices)
        
        self.assertFloatEqual(test_stat_result, test_stat_output)
        self.assertFloatEqual(pvals_result, pvals_output)
        self.assertEqual(means_result, means_output)

    def test_fdr_correction(self):
        """fdr_correction works"""
        pvals_1 = [0.56370286165077377,
             0.043308142810792101,
             0.38363032713198986,
             0.08326451666355042,
             0.043308142810792101,
             0.77282999268444919]
        pvals_1_output = [0.67644343398092854,
             0.25984885686475262,
             0.57544549069798479,
             0.16652903332710084,
             0.12992442843237631,
             0.77282999268444919]

        pvals_1_result = fdr_correction(pvals_1)
        self.assertFloatEqual(pvals_1_result, pvals_1_output)

    def test_bonferroni_correction(self):
        """bonferroni_correction works"""
        pvals_1 = [0.56370286165077377,
             0.043308142810792101,
             0.38363032713198986,
             0.08326451666355042,
             0.043308142810792101,
             0.77282999268444919]
        pvals_1_output = [3.3822171699046426,
             0.2598488568647526,
             2.301781962791939,
             0.4995870999813025,
             0.2598488568647526,
             4.636979956106695]
        
        pvals_1_result = bonferroni_correction(pvals_1)
        self.assertFloatEqual(pvals_1_result, pvals_1_output)  

    def test_output_formatter(self):
        """output_formatter works"""
        #Using ANOVA test for example
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""
        bt = parse_biom_table_str(otu_table_1)

        test_stats = [0.18990151199889027,
             6.7142857142857144,
             1.8424031345232912,
             1.4643841007477372,
             4.5675332910589734,
             9.8389688760617899e-05]
        pvals = [0.6782397284636259,
             0.041145883121579234,
             0.22350244183135481,
             0.27174025956151771,
             0.076447615888438403,
             0.9924073718332751]
        means = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]
        fdr_pvals = [0.8138876741563511,
             0.24687529872947539,
             0.44700488366270963,
             0.40761038934227656,
             0.22934284766531521,
             0.9924073718332751]
        bon_pvals = [4.069438370781755,
             0.2468752987294754,
             1.3410146509881289,
             1.6304415573691062,
             0.4586856953306304,
             5.95444423099965]
        cat_sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}

        line_header_out = 'OTU\tTest-Statistic\tP\tFDR_P\tBonferroni_P\tcat1_mean\tcat2_mean\tTaxonomy'
        line_1_out = 'OTU1\t0.189901511999\t0.678239728464\t0.813887674156\t4.06943837078\t52.25\t43.0\tk__One'
        line_6_out = 'OTU6\t9.83896887606e-05\t0.992407371833\t0.992407371833\t5.954444231\t48.0\t47.75\tk__Six'

        lines = output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals, \
            means, cat_sample_indices)

        self.assertEqual(lines[0], line_header_out)
        self.assertEqual(lines[1], line_1_out)
        self.assertEqual(lines[6], line_6_out)






#run unit tests if run from command-line
if __name__ == '__main__':
    main()













    
