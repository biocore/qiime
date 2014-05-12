#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"

import tempfile
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.util import get_qiime_temp_dir
from qiime.split_libraries_lea_seq import (get_cluster_ratio, get_consensus,
                                         read_input_file, extract_primer,
                                         select_unique_rand_bcs)
import os


class WorkflowTests(TestCase):

    def setUp(self):
        """setup the test values"""
        temp_dir = get_qiime_temp_dir()
        self.fasta_seqs_of_rand_bcs = fasta_seqs_of_rand_bcs
        self.fasta_seqs_for_uclust = fasta_seqs_for_uclust
        self.fasta_seqs_for_consensus = fasta_seqs_for_consensus
        self.fasta_file_for_consensus = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=temp_dir)
        self.fasta_file_for_cluster_ratio = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=temp_dir)

    def tearDown(self):
        """remove all the files after completing tests """
        os.unlink(self.fasta_file_for_consensus.name)
        os.unlink(self.fasta_file_for_cluster_ratio.name)
        # delete = False in creation of tempfiles
        # as they need to be opened in functions
        # hence the need to remove them

    def test_select_unique_rand_bcs(self):
        fasta_seqs_of_rand_bcs = self.fasta_seqs_of_rand_bcs
        actual = select_unique_rand_bcs(fasta_seqs_of_rand_bcs)
        expected = {'ATTGCATTGCATTGCATTGC': 1, 'ATTGCATTGCATTGCATTG': 0, 'ATTGCTTATTGCATTGCTTT': 1}
        self.assertEqual(actual, expected)

    def test_get_consensus(self):
        fasta_seqs_for_consensus = self.fasta_seqs_for_consensus
        temp_file = self.fasta_file_for_consensus
        temp_file.write(fasta_seqs_for_consensus)
        temp_file_name = temp_file.name
        temp_file.close()
        temp_file = open(temp_file_name, 'r')
        actual = get_consensus(temp_file, 2)
        expected = 'ATGCATGC'
        self.assertEqual(actual, expected)

    def test_get_cluster_ratio(self):
        fasta_seqs_for_uclust = self.fasta_seqs_for_uclust
        temp_dir = get_qiime_temp_dir()
        temp_file = self.fasta_file_for_cluster_ratio
        temp_file.write(fasta_seqs_for_uclust)
        temp_file_name = temp_file.name
        temp_file.close()
        actual = get_cluster_ratio(temp_file_name)
        expected = 0.125
        self.assertEqual(actual, expected)

    def test_read_input_file(self):
        pass

    def test_extract_primers(self):
        pass

fasta_seqs_for_uclust = """>abc|1
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
"""

fasta_seqs_for_consensus=""">id1|1
ATGCATGG
>id2|14
ATGCATGC
"""

fasta_seqs_of_rand_bcs = ['ATTGCATTGCATTGCATTGC', 'ATTGCATTGCATTGCATTGC', 'ATTGCATTGCATTGCATTG', 'ATTGCTTATTGCATTGCTTT']

# run tests if called from command line
if __name__ == '__main__':
    main()
