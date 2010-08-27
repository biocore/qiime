#!/usr/bin/env python
# File created on 17 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"
 
from os.path import exists
from os import rmdir
from numpy import array
from matplotlib.axes import Subplot

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from cogent.app.util import get_tmp_filename

from qiime.plot_rank_abundance_graph import make_sorted_frequencies,\
    plot_rank_abundance_graph, plot_rank_abundance_graphs
from qiime.util import create_dir

class PlotRankAbundance(TestCase):
    
    def setUp(self):

        self.tmp_dir = None
        self.files_to_remove = []
    
    def tearDown(self):
        
        remove_files(self.files_to_remove)
        if self.tmp_dir:
            rmdir(self.tmp_dir)
        
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
 
        self.otu_table = otu_table_fake.split('\n')       
        self.tmp_dir = get_tmp_filename(prefix="test_plot_rank_abundance",
                                   suffix="/")
        create_dir(self.tmp_dir)


    
        #test all supported filetypes
        for file_type in ['pdf','svg','png','eps']:
            plot_rank_abundance_graphs('S3', iter(self.otu_table), self.tmp_dir,
                                       file_type=file_type)
            tmp_file = self.tmp_dir+"rank_abundanceS3."+file_type
            self.files_to_remove.append(tmp_file)
            self.assertTrue(exists(tmp_file))
            
    def test_plot_rank_abundance_graphs(self):
        """plot_rank_abundance_graphs works with any number of samples"""
 
        self.otu_table = otu_table_fake.split('\n')       
        self.tmp_dir = get_tmp_filename(prefix="test_plot_rank_abundance",
                                   suffix="/")
        create_dir(self.tmp_dir)
 
        #test empty sample name
        self.assertRaises(ValueError, plot_rank_abundance_graphs, '',
                          iter(self.otu_table), self.tmp_dir)
        #test invalid sample name
        self.assertRaises(ValueError, plot_rank_abundance_graphs,
                          'Invalid_sample_name',
                          iter(self.otu_table), self.tmp_dir)

        #test with two samples
        file_type="pdf"
        plot_rank_abundance_graphs('S3,S5', iter(self.otu_table), self.tmp_dir,
                                       file_type=file_type)
        tmp_file = self.tmp_dir+"rank_abundanceS3_S5."+file_type
        self.files_to_remove.append(tmp_file)
        self.assertTrue(exists(tmp_file)) 
 
        # test with all samples
        plot_rank_abundance_graphs('*', iter(self.otu_table), self.tmp_dir,
                                       file_type=file_type)
        tmp_file = self.tmp_dir+"rank_abundanceS3_S4_S5."+file_type
        self.files_to_remove.append(tmp_file)
        self.assertTrue(exists(tmp_file)) 


otu_table_fake = """#Full OTU Counts
#OTU ID	S3	S4	S5	Consensus Lineage
0	1	0	1	Root;Bacteria
3	2	0	1	Root;Bacteria;Acidobacteria
4	1	0	9	Root;Bacteria;Bacteroidetes
2	1	0	1	Root;Bacteria;Acidobacteria;Acidobacteria;Gp5
6	1	25	42	Root;Archaea"""  

if __name__ == "__main__":
    main()
