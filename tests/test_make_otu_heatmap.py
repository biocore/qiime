#!/usr/bin/env python
#file test_make_otu_heatmap_html.py

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Release"


from numpy import array, log10, asarray, float64
from cogent.util.unit_test import TestCase, main
from qiime.make_otu_heatmap import (extract_metadata_column, 
    get_order_from_categories, get_order_from_tree, make_otu_labels, 
    names_to_indices, get_log_transform, get_clusters, 
    get_fontsize, plot_heatmap)
from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files
from os.path import exists

class TopLevelTests(TestCase):
    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""
        self.col_header=['Sample1', 'Sample2', 'Sample3', 
                         'Sample4', 'Sample5', 'Sample6']
        self.row_header=['OTU1','OTU2','OTU3']
        self.otu_table=array([[0,0,9,5,3,1],
                              [1,5,4,0,3,2],
                              [2,3,1,1,2,5]])
        self.lineages=[['Bacteria'],['Archaea'],['Streptococcus']]
        self.full_lineages=[['1A','1B','1C','Bacteria'],
                            ['2A','2B','2C','Archaea'],
                            ['3A','3B','3C','Streptococcus']]
        self.metadata = [[['Sample1','NA','A'],
                          ['Sample2','NA','B'],
                          ['Sample3','NA','A'],
                          ['Sample4','NA','B'],
                          ['Sample5','NA','A'],
                          ['Sample6','NA','B']],
                          ['SampleID', 'CAT1','CAT2'],[]]
        self.tree_text = ["('OTU3',('OTU1','OTU2'))"]
        self.tmp_heatmap_fpath = get_tmp_filename(
            prefix='test_heatmap_',
            suffix='.pdf'
            )


    def test_extract_metadata_column(self):
        """Extracts correct column from mapping file"""
        obs = extract_metadata_column(self.col_header, 
                self.metadata, category='CAT2')
        exp = ['A','B','A','B','A','B']
        self.assertEqual(obs,exp)
        
    def test_get_order_from_categories(self):
        """Sample indices should be clustered within each category"""
        category_labels = ['A','B','A','B','A','B']
        obs = get_order_from_categories(self.otu_table, category_labels)
        exp = [0,4,2,1,5,3]
        self.assertEqual(obs,exp)

    def test_get_order_from_tree(self):
        obs = get_order_from_tree(self.row_header, self.tree_text)
        exp = [2,0,1]
        self.assertEqual(obs,exp)
        
    def test_make_otu_labels(self):
        obs = make_otu_labels(self.row_header, self.lineages, n_levels=1)
        exp = ['Bacteria (OTU1)', 'Archaea (OTU2)', 'Streptococcus (OTU3)']
        self.assertEqual(obs,exp)

        obs = make_otu_labels(self.row_header, self.full_lineages, n_levels=3)
        exp = ['1B;1C;Bacteria (OTU1)', 
               '2B;2C;Archaea (OTU2)',
               '3B;3C;Streptococcus (OTU3)']
        self.assertEqual(obs,exp)
            
    def test_names_to_indices(self):
        new_order = ['Sample4', 'Sample2', 'Sample3', 
                     'Sample6', 'Sample5', 'Sample1']
        obs = names_to_indices(self.col_header, new_order)
        exp = [3, 1, 2, 5, 4, 0]
        self.assertEqual(obs,exp)

    def test_get_log_transform(self):
        eps = .01
        obs = get_log_transform(self.otu_table, eps=eps)
        xform = asarray(self.otu_table, dtype=float64)
        xform[xform==0] = eps
    
        for i, row in enumerate(obs):
            self.assertEqual(row, log10(xform[i]))

    def test_get_clusters(self):
        obs = get_clusters(self.otu_table, axis='row')
        exp = [0, 1, 2]
        self.assertEqual(obs,exp)
        obs = get_clusters(self.otu_table, axis='column')
        exp = [0, 5, 4, 1, 2, 3]
        self.assertEqual(obs,exp)
        
    def test_plot_heatmap(self):
        plot_heatmap(self.otu_table, self.row_header, self.col_header, 
                     filename=self.tmp_heatmap_fpath)
        self.assertEqual(exists(self.tmp_heatmap_fpath), True)
        remove_files(set([self.tmp_heatmap_fpath]))


#run tests if called from command line
if __name__ == "__main__":
    main()
