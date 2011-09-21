#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@gmail.com"
__status__ = "Development"

from cogent.util.misc import remove_files

from random import seed
from unittest import TestCase, main

from qiime.subsample_fasta import subsample_fasta
from qiime.util import get_tmp_filename, load_qiime_config
    
class SubSampleFastaTests(TestCase):
    """ """
    
    def setUp(self):
        """ """
        self.expected_lines_50_perc = expected_lines_50_perc
        self.expected_lines_20_perc = expected_lines_20_perc
        
        self.temp_dir = load_qiime_config()['temp_dir']
        
        self.fasta_lines = fasta_lines
        self.fasta_filepath = get_tmp_filename(
            prefix='subsample_test_', suffix='.fasta')
        self.fasta_file = open(self.fasta_filepath, "w")
        self.fasta_file.write(self.fasta_lines)
        self.fasta_file.close()
        
        self.output_filepath = get_tmp_filename(prefix='subsample_output_',
         suffix='.fasta')
        
        self._files_to_remove =\
         [self.fasta_filepath]

    def tearDown(self):
        remove_files(self._files_to_remove)
        
    def test_subsample_fasta_50(self):
        """ subsample_fasta correctly subsamples input fasta file """
        
        # fixed seed for consistent calls with random()
        seed(128)
        
        subsample_fasta(self.fasta_filepath, self.output_filepath,
         percent_subsample = 0.50)
    
        self._files_to_remove.append(self.output_filepath)
         
        actual_results =\
         [line.strip() for line in open(self.output_filepath, "U")]
        
        self.assertEqual(actual_results, self.expected_lines_50_perc)

    def test_subsample_fasta_20(self):
        """ subsample_fasta correctly subsamples input fasta file """
        
        seed(12210)

        subsample_fasta(self.fasta_filepath, self.output_filepath,
         percent_subsample = 0.20)
    
        self._files_to_remove.append(self.output_filepath)
         
        actual_results =\
         [line.strip() for line in open(self.output_filepath, "U")]
        
        self.assertEqual(actual_results, self.expected_lines_20_perc)
        
        
        
        
    


# Long strings of test data go here
fasta_lines = """>seq1
ACCAGCGGAGAC
>seq2
ACAGAGAGACCC
>seq3
ATTACCAGATTAC
>seq4
ACAGGAGACCGAGAAGA
>seq5
ACCAGAGACCGAGA
"""

expected_lines_50_perc = """>seq2
ACAGAGAGACCC
>seq4
ACAGGAGACCGAGAAGA
>seq5
ACCAGAGACCGAGA""".split('\n')

expected_lines_20_perc = """>seq4
ACAGGAGACCGAGAAGA""".split('\n')

if __name__ == "__main__":
    main()
