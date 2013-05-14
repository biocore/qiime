#!/usr/bin/env python
# File created on 17 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder","Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"
 
from os.path import exists,abspath
from shutil import rmtree
from numpy import array
from matplotlib.axes import Subplot

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.util import get_tmp_filename
from qiime.plot_rank_abundance_graph import make_sorted_frequencies,\
    plot_rank_abundance_graph, plot_rank_abundance_graphs
from qiime.util import create_dir
from biom.parse import parse_biom_table_str

class PlotRankAbundance(TestCase):
    
    def setUp(self):
        self.tmp_dir='/tmp/'
        self.files_to_remove = []
        self._dirs_to_remove=[]
        
    def tearDown(self):
        
        remove_files(self.files_to_remove)
        if self._dirs_to_remove:
            for i in self._dirs_to_remove:
                rmtree(i)
        
    def test_make_sorted_frequencies(self):
        """make_sorted_frequencies transforms and sorts correctly"""
        
        #works on empty
        counts = array([])
        self.assertEqual(make_sorted_frequencies(counts),[])

         #works on zeros
        counts = array([0,0,0,0,0,0])
        self.assertEqual(make_sorted_frequencies(counts),[])

        #works on flat data
        counts = array([3,3,3,3,3])
        expected_freqs = [0.2, 0.2, 0.2, 0.2, 0.2] 
        observed_freqs = make_sorted_frequencies(counts)
        self.assertEqual(observed_freqs, expected_freqs)

        #works on real data
        counts = array([1,2,0,1,0,2,4])
        expected_freqs = [0.4, 0.2, 0.2, 0.1, 0.1] 
        observed_freqs = make_sorted_frequencies(counts)
        self.assertEqual(observed_freqs, expected_freqs)
   
    def test_make_sorted_frequencies_abolute(self):
        """make_sorted_frequencies returns correct absolute values"""
        
        #works on empty
        counts = array([])
        self.assertEqual(make_sorted_frequencies(counts,True),[])
  
        #works on zeros
        counts = array([0,0,0,0,0,0])
        self.assertEqual(make_sorted_frequencies(counts,True),[])
        
        #works on flat data
        counts = array([3,3,3,3,3])
        expected_freqs = [3,3,3,3,3] 
        observed_freqs = make_sorted_frequencies(counts,True)
        self.assertEqual(observed_freqs, expected_freqs)

        #works o real data
        counts = array([1,2,0,1,0,2,4])
        expected_freqs = [4,2,2,1,1] 
        observed_freqs = make_sorted_frequencies(counts,True)
        self.assertEqual(observed_freqs, expected_freqs)

    def test_plot_rank_abundance_graph(self):
        """plot_rank_abudance_graph plots something"""

        counts = array([20,15,12,8,4,2,1,3,1,2])
        observed = plot_rank_abundance_graph(counts)
        
        # can we test something more clever here?
        # basically just tests that something is drawn, but not what
        self.assertEqual(type(observed), Subplot)


    def test_plot_rank_abundance_graphs_filetype(self):
        """plot_rank_abundance_graphs works with all filetypes"""
 
        self.otu_table = parse_biom_table_str(otu_table_sparse)
        self.dir = get_tmp_filename(tmp_dir=self.tmp_dir,
                                   prefix="test_plot_rank_abundance",
                                   suffix="/")
                                   
        create_dir(self.dir)
        self._dirs_to_remove.append(self.dir)
    
        #test all supported filetypes
        for file_type in ['pdf','svg','png','eps']:
            tmp_file = abspath(self.dir+"rank_abundance_cols_0."+file_type)

            plot_rank_abundance_graphs(tmp_file,'S3', self.otu_table, file_type=file_type)
            self.files_to_remove.append(tmp_file)
            self.assertTrue(exists(tmp_file))
            
    def test_plot_rank_abundance_graphs_sparse(self):
        """plot_rank_abundance_graphs works with any number of samples (SparseOTUTable)"""
 
        self.otu_table = parse_biom_table_str(otu_table_sparse)
        self.dir = get_tmp_filename(tmp_dir=self.tmp_dir,
                                   prefix="test_plot_rank_abundance",
                                   suffix="/")
        create_dir(self.dir)
        self._dirs_to_remove.append(self.dir)
        tmp_fname = get_tmp_filename(tmp_dir=self.dir)

        #test empty sample name
        self.assertRaises(ValueError, plot_rank_abundance_graphs, tmp_fname, '',
                          self.otu_table)
        #test invalid sample name
        self.assertRaises(ValueError, plot_rank_abundance_graphs, tmp_fname,
                          'Invalid_sample_name',
                          self.otu_table)

        #test with two samples
        file_type="pdf"
        tmp_file = abspath(self.dir+"rank_abundance_cols_0_2."+file_type)

        plot_rank_abundance_graphs(tmp_file, 'S3,S5', self.otu_table,
                                       file_type=file_type)

        self.assertTrue(exists(tmp_file)) 
        self.files_to_remove.append(tmp_file)
        # test with all samples
        tmp_file = abspath(self.dir+"rank_abundance_cols_0_1_2."+file_type)
        plot_rank_abundance_graphs(tmp_file, '*', self.otu_table,
                                       file_type=file_type)
        
        self.files_to_remove.append(tmp_file)
        self.assertTrue(exists(tmp_file)) 

    def test_plot_rank_abundance_graphs_dense(self):
        """plot_rank_abundance_graphs works with any number of samples (DenseOTUTable)"""
 
        self.otu_table = parse_biom_table_str(otu_table_dense)
        self.dir = get_tmp_filename(tmp_dir=self.tmp_dir,
                                   prefix="test_plot_rank_abundance",
                                   suffix="/")
        create_dir(self.dir)
        self._dirs_to_remove.append(self.dir)
        tmp_fname = get_tmp_filename(tmp_dir=self.dir)

        #test empty sample name
        self.assertRaises(ValueError, plot_rank_abundance_graphs, tmp_fname,'',
                          self.otu_table)
        #test invalid sample name
        self.assertRaises(ValueError, plot_rank_abundance_graphs, tmp_fname,
                          'Invalid_sample_name',
                          self.otu_table)

        #test with two samples
        file_type="pdf"
        tmp_file = abspath(self.dir+"rank_abundance_cols_0_2."+file_type)
        plot_rank_abundance_graphs(tmp_file, 'S3,S5', self.otu_table,
                                       file_type=file_type)

        self.assertTrue(exists(tmp_file)) 
        self.files_to_remove.append(tmp_file)
        # test with all samples
        tmp_file = abspath(self.dir+"rank_abundance_cols_0_1_2."+file_type)

        plot_rank_abundance_graphs(tmp_file,'*', self.otu_table,file_type=file_type)
        
        self.files_to_remove.append(tmp_file)
        self.assertTrue(exists(tmp_file)) 


otu_table_sparse = """{"rows": [{"id": "0", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "3", "metadata": {"taxonomy": ["Root", "Bacteria", "Acidobacteria"]}}, {"id": "4", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "2", "metadata": {"taxonomy": ["Root", "Bacteria", "Acidobacteria", "Acidobacteria", "Gp5"]}}, {"id": "6", "metadata": {"taxonomy": ["Root", "Archaea"]}}], "format": "Biological Observation Matrix v0.9", "data": [[0, 0, 1.0], [0, 2, 1.0], [1, 0, 2.0], [1, 2, 1.0], [2, 0, 1.0], [2, 2, 9.0], [3, 0, 1.0], [3, 2, 1.0], [4, 0, 1.0], [4, 1, 25.0], [4, 2, 42.0]], "columns": [{"id": "S3", "metadata": null}, {"id": "S4", "metadata": null}, {"id": "S5", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2571", "matrix_type": "sparse", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T19:33:37.780300", "type": "OTU table", "id": null, "matrix_element_type": "float"}"""

otu_table_dense = """{"rows": [{"id": "0", "metadata": {"taxonomy": ["Root", "Bacteria"]}}, {"id": "3", "metadata": {"taxonomy": ["Root", "Bacteria", "Acidobacteria"]}}, {"id": "4", "metadata": {"taxonomy": ["Root", "Bacteria", "Bacteroidetes"]}}, {"id": "2", "metadata": {"taxonomy": ["Root", "Bacteria", "Acidobacteria", "Acidobacteria", "Gp5"]}}, {"id": "6", "metadata": {"taxonomy": ["Root", "Archaea"]}}], "format": "Biological Observation Matrix v0.9", "data": [[1, 0, 1], [2, 0, 1], [1, 0, 9], [1, 0, 1], [1, 25, 42]], "columns": [{"id": "S3", "metadata": null}, {"id": "S4", "metadata": null}, {"id": "S5", "metadata": null}], "generated_by": "QIIME 1.4.0-dev, svn revision 2571", "matrix_type": "dense", "shape": [5, 3], "format_url": "http://www.qiime.org/svn_documentation/documentation/biom_format.html", "date": "2011-12-21T19:33:28.922480", "type": "OTU table", "id": null, "matrix_element_type": "int"}"""

if __name__ == "__main__":
    main()
