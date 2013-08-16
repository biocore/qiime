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
from numpy import array
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

    def test_run_ocs_test(self):
        """run_ocs_test works"""
        




#run unit tests if run from command-line
if __name__ == '__main__':
    main()













    
