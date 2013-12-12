#!/usr/bin/env python

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Development"

from os.path import isdir, isfile, exists
from shutil import rmtree

from cogent.util.unit_test import TestCase, main
from qiime.util import get_tmp_filename
from cogent.util.misc import remove_files, get_random_directory_name

from qiime.quality_scores_plot import generate_histogram,\
 plot_qual_report, get_qual_stats, bin_qual_scores, write_qual_report
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
        
        self.output_dir = get_random_directory_name(prefix = '/tmp/')
        self.output_dir += '/'
        
        create_dir(self.output_dir)
        
        self.expected_output_text_file = expected_output_text_file
        
        
        self._files_to_remove =\
         [self.qual_fp]
        
    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir ):
            rmtree(self.output_dir )
            
    def test_generate_histogram(self):
        """ No errors when calling function, creates output files"""
        
        # Cannot test content of graphics file, only successful execution
        
        output_dir = self.output_dir
        
        
        # Should not raise an error with good data
        
        generate_histogram(self.qual_fp, output_dir, verbose = False)
        
        expected_outfile = output_dir + 'quality_scores_plot.pdf'
        
        self.assertTrue(isfile(expected_outfile))
        
        # Test text file output for proper data
        text_output_fp = output_dir + "quality_bins.txt"
        
        text_output_f = open(text_output_fp, "U")
        
        actual_text_output = "\n".join([line.strip() for line in text_output_f])


        self.assertEqual(actual_text_output, self.expected_output_text_file)
        
    def test_write_qual_report(self):
        """ Writes data to output text file properly """
        
        
        output_dir = self.output_dir
        
        
        qual_bins = [[1, 2, 6], [1, 2, 3], [2, 4], [4]]
        
        expected_ave_bins = [3.00, 2.00, 3.00, 4.00]
        expected_std_dev_bins = [2.16, 0.816, 1.0, 0]
        expected_total_bases_bins = [3, 3, 2, 1]
        
        score_min = 25
        
        write_qual_report(expected_ave_bins, expected_std_dev_bins,
         expected_total_bases_bins, output_dir, score_min)
        
        # Test text file output for proper data
        text_output_fp = output_dir + "quality_bins.txt"
        
        text_output_f = open(text_output_fp, "U")
        
        actual_text_output = [line.strip() for line in text_output_f]
        
        ave_bin_index = 2
        std_dev_bin_index = 4
        total_bases_index = 6
        
        actual_bins_ave =\
         [float(f) for f in actual_text_output[ave_bin_index].split(',')]
        actual_bins_std =\
         [float(f) for f in actual_text_output[std_dev_bin_index].split(',')]
        actual_bins_total_bases =\
         [float(f) for f in actual_text_output[total_bases_index].split(',')]
        

        self.assertEqual(actual_bins_ave, expected_ave_bins)
        self.assertEqual(actual_bins_std, expected_std_dev_bins)
        self.assertEqual(actual_bins_total_bases, expected_total_bases_bins)
        
        
        
        
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
        score_min = 25
        
        actual_ave_bins, actual_std_dev_bins, actual_total_bases_bins,\
         suggested_trunc_pos = get_qual_stats(qual_bins, score_min)
         
        # Should give correct suggested truncation position, where the quality
        # score average went below 25, in this case, at base 0
        expected_trunc_pos = 0
        
        self.assertEqual(suggested_trunc_pos, expected_trunc_pos)
         
        # Round standard deviation calculations
        for n in range(len(actual_std_dev_bins)):
            actual_std_dev_bins[n] = round(actual_std_dev_bins[n], 3)
         
        self.assertEqual(actual_ave_bins, expected_ave_bins)
        self.assertEqual(actual_std_dev_bins, expected_std_dev_bins)
        self.assertEqual(actual_total_bases_bins, expected_total_bases_bins)
        
        
    def test_plot_qual_report(self):
        """ Is called without error """
        
        output_dir = self.output_dir
        
        
        score_min = 20
        
        ave_bins = [3, 2, 3, 4]
        std_dev_bins = [2.16, 0.816, 1.0, 0]
        total_bases_bins = [3, 3, 2, 1]
        
        plot_qual_report(ave_bins, std_dev_bins, total_bases_bins, score_min,
         output_dir)
        
        
        
        
        
        
        
        
        
        
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

expected_output_text_file = """# Suggested nucleotide truncation position (None if quality score average did not drop below the score minimum threshold): 245
# Average quality score bins
34.333,34.333,34.333,34.667,35.000,34.333,34.333,34.333,35.000,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,37.667,37.667,38.000,38.000,31.667,30.667,30.667,35.333,35.333,36.000,35.333,35.333,36.000,31.333,26.000,26.000,28.667,31.667,30.667,28.333,32.000,32.000,35.000,36.667,36.667,37.333,38.000,35.667,35.000,35.000,35.000,35.667,33.000,33.000,33.000,35.333,35.333,35.667,35.667,35.667,34.667,34.667,33.333,33.333,33.333,35.000,34.667,34.667,34.667,34.667,35.667,35.667,35.667,35.667,35.667,35.667,34.333,34.333,34.333,32.667,31.667,31.667,31.667,34.667,35.667,35.667,35.667,35.667,35.667,35.667,35.000,34.333,34.333,34.333,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,34.333,34.333,34.000,35.333,32.667,32.667,32.000,32.000,34.000,35.000,35.000,35.667,35.667,35.667,35.667,35.333,35.333,35.333,35.333,35.333,35.333,35.333,35.667,35.667,35.667,35.667,35.667,35.667,34.000,32.667,32.667,29.667,31.000,32.000,32.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.000,34.667,34.667,35.333,34.667,33.000,33.000,32.333,32.333,32.333,33.000,28.667,28.667,28.667,28.667,28.667,31.000,31.000,34.000,34.000,35.667,34.000,31.333,31.000,30.333,30.333,30.000,34.000,35.667,34.667,32.333,32.333,29.000,30.667,30.333,30.333,30.333,32.667,32.000,32.000,33.667,33.667,34.333,35.333,35.333,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,35.667,33.667,32.000,32.000,32.000,33.667,35.667,33.667,33.667,32.667,34.667,35.000,31.667,31.333,32.333,32.333,35.000,35.000,35.000,33.333,32.333,27.333,27.333,26.667,26.667,26.667,30.000,29.000,32.333,32.333,34.000,34.000,34.333,34.000,34.000,34.667,35.000,32.000,30.667,29.333,29.333,24.667,28.000,28.000,28.000,26.000,25.333,27.667,29.000,29.000,29.000,30.000,30.000,29.500,29.500,28.500,26.000,23.500,23.500,26.500,28.500,29.500,25.500,24.000,23.500,22.000,30.000,30.000
# Standard deviation bins
1.700,1.700,1.700,2.055,1.633,2.494,2.494,2.494,1.633,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,2.055,2.055,2.160,2.160,6.944,8.219,8.219,3.682,3.682,2.944,2.055,2.055,2.160,6.944,10.033,10.033,8.179,8.498,8.219,8.498,6.164,6.164,3.559,2.357,2.357,2.055,2.160,0.943,1.633,1.633,1.633,0.943,1.633,1.633,1.633,0.471,0.471,0.943,0.943,0.943,2.055,2.055,2.625,2.625,2.625,1.633,2.055,2.055,2.055,2.055,0.943,0.943,0.943,0.943,0.943,0.943,2.494,2.494,2.494,4.784,4.497,4.497,4.497,2.055,0.943,0.943,0.943,0.943,0.943,0.943,1.633,2.494,2.494,2.494,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,2.494,2.494,2.449,1.247,4.784,4.784,4.546,4.546,2.160,1.633,1.633,0.943,0.943,0.943,0.943,1.247,1.247,1.247,1.247,1.247,1.247,1.247,0.943,0.943,0.943,0.943,0.943,0.943,2.944,4.784,4.784,5.185,4.546,4.546,4.784,0.943,0.943,0.943,0.943,0.943,0.943,0.943,1.633,1.700,1.700,1.247,1.700,4.320,4.320,5.249,5.249,5.249,4.320,10.403,10.403,10.403,10.403,10.403,4.320,4.320,0.816,0.816,0.943,2.944,6.650,5.354,6.799,6.799,6.683,2.944,0.943,2.055,5.249,5.249,9.933,7.587,7.318,7.318,3.300,3.300,2.944,2.944,0.943,0.943,0.943,0.471,0.471,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,0.943,3.399,5.715,5.715,5.715,3.399,0.943,3.399,3.399,4.784,1.700,1.633,4.110,3.682,3.859,3.859,0.816,0.816,0.816,1.700,1.886,6.128,6.128,4.989,4.989,4.989,2.160,5.715,1.886,1.886,2.160,2.160,1.700,2.160,2.160,1.247,0.816,4.967,4.497,6.799,6.799,5.312,2.449,2.449,2.449,3.742,4.497,4.028,1.000,1.000,1.000,0.000,0.000,0.500,0.500,0.500,4.000,1.500,1.500,1.500,0.500,0.500,0.500,2.000,1.500,0.000,0.000,0.000
# Total bases per nucleotide position bins
3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1"""

if __name__ == "__main__":
    main()
