#!/usr/bin/env python

__author__ = "William Walters"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Development"

from os.path import isdir, isfile
from shutil import rmtree

from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
from cogent.util.misc import remove_files

from qiime.quality_scores_plot import generate_histogram,\
 plot_qual_report, get_qual_stats, bin_qual_scores
from qiime.util import create_dir

class QualityScoresPlotTests(TestCase):
    """ Unit tests for the quality_scores_histogram.py module """
    
    def setUp(self):
        # create the temporary input files that will be used
        
        self._files_to_remove = []
        
        
        self.qual_fp = get_tmp_filename(\
         prefix = 'qual_scores_',
         suffix = '.qual')
        seq_file = open(self.qual_fp, 'w')
        seq_file.write(qual_scores)
        seq_file.close()
        
        
        self._files_to_remove =\
         [self.qual_fp]
        
    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if isdir('/tmp/test_dir_qual_scores/'):
            rmtree('/tmp/test_dir_qual_scores/')
            
    def test_generate_histogram(self):
        """ No errors when calling function, creates output file"""
        
        # Cannot test content of graphics file, only successful execution
        
        output_dir = '/tmp/test_dir_qual_scores/'
        
        create_dir(output_dir)
        
        # Should not raise an error with good data
        
        generate_histogram(self.qual_fp, output_dir)
        
        
        
    def test_bin_qual_scores(self):
        """ Properly bins qual scores according to nucleotide position """
        
        qual_data = {'seq1':[10, 20, 30, 40], 'seq2':[11, 21, 31],
         'seq3':[12, 22]}
        
         
        expected_bins = [[10, 11, 12], [20, 21, 22], [30, 31], [40]]
        
        actual_bins = bin_qual_scores(qual_data)
        
        
        # Because of arbritrary dictionary order, need to sort results
        
        for bin in range(len(actual_bins)):
            actual_bins[bin].sort()
        
        self.assertEqual(actual_bins, expected_bins)
        
        
    def test_get_qual_stats(self):
        """ Properly generates averages, std dev for bins """
        
        qual_bins = [[1, 2, 6], [1, 2, 3], [2, 4], [4]]
        
        expected_ave_bins = [3, 2, 3, 4]
        expected_std_dev_bins = [2.16, 0.816, 1.0, 0]
        expected_total_bases_bins = [3, 3, 2, 1]
        
        actual_ave_bins, actual_std_dev_bins, actual_total_bases_bins =\
         get_qual_stats(qual_bins)
         
        # Round standard deviation calculations
        for n in range(len(actual_std_dev_bins)):
            actual_std_dev_bins[n] = round(actual_std_dev_bins[n], 3)
         
        self.assertEqual(actual_ave_bins, expected_ave_bins)
        self.assertEqual(actual_std_dev_bins, expected_std_dev_bins)
        self.assertEqual(actual_total_bases_bins, expected_total_bases_bins)
        
        
    def test_plot_qual_report(self):
        """ Is called without error, creates output file """
        
        output_dir = '/tmp/test_dir_qual_scores/'
        
        create_dir(output_dir)
        
        score_min = 20
        
        ave_bins = [3, 2, 3, 4]
        std_dev_bins = [2.16, 0.816, 1.0, 0]
        total_bases_bins = [3, 3, 2, 1]
        
        plot_qual_report(ave_bins, std_dev_bins, total_bases_bins, score_min,
         output_dir)
        
        expected_outfile = '/tmp/test_dir_qual_scores/quality_scores_plot.pdf'
        
        self.assertTrue(isfile(expected_outfile))
        
        
        
        
        
        
# Long strings at the end for better readability
   


qual_scores = """>seq1
36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 40 40 40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37
37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36
36 36 36 36 36 31 31 31 31 31 31 31
>seq2
35 35 35 35 33 31 31 31 33 35 35 35 35 35 35 35 35 35 35 35 35 35 23 20 20 31 31 33 33 33 35 23 17 17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35
31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35
35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27 14 14 14 14 14 25 25 34 34 35 35 35 32 33 33 32 35 35 32 25 25
15 20 20 20 28 35 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 29 24 24 24 29 35 35 35 35 33 33 31 31 34 34 34 34 34 34 31 20 20 20 20 20 31 34 31 31 31 31 32 31 31 33
34 25 25 20 20 18 25 28 28 22 20 22 28 28 28 30 30 29 29 29 30 25 25 25 29 29 26 26 25
>seq3
32 32 32 32 35 35 35 35 35 35 35 35 35 35 35 35 35 35 38 38 39 39 32 32 32 35 35 35 35 35 34 31 21 21 25 35 32 25 25 25 32 35 35 37 39 35 35 35 35 35 35 35 35 35 35 35 35 35 32 32
32 32 32 35 32 32 32 32 35 35 35 35 35 35 35 35 35 35 32 32 32 32 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 34 34 26 26 26 26 32 35 35 35 35 35
35 34 34 34 34 34 34 34 35 35 35 35 35 35 35 35 35 26 26 26 26 35 35 35 35 35 35 35 35 34 34 34 34 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 22 24 21 21 21 30 35 35 35 35
35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 29 29 26 34 35 27 27 27 27 35 35 35 35 35 35 35 32 32 32 32 32 35 35 35 35 35 35 35 35
35 35 31 32 32 25 28 25 25 25 25 30 30 30 30 30 30 30 30 28 22 22 22 28 28 30 25 22 22 22 30 30"""


if __name__ == "__main__":
    main()
