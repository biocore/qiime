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
                                         get_LEA_seq_consensus_seqs, extract_primer,
                                         select_unique_rand_bcs, SeqLengthMismatchError)
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
        self.fasta_file_for_consensus_unequal_length = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=temp_dir)
        self.fasta_file_for_consensus_name = self.fasta_file_for_consensus.name
        self.fasta_file_for_cluster_ratio_name = self.fasta_file_for_cluster_ratio.name
        self.fasta_file_for_consensus_unequal_length_name = self.fasta_file_for_consensus_unequal_length.name
        self.fasta_seq_for_primer = 'AATGCCCCC'
        self.possible_primers = ['ATGC', 'ATTT']


    def tearDown(self):
        """remove all the files after completing tests """
        os.unlink(self.fasta_file_for_consensus_name)
        os.unlink(self.fasta_file_for_cluster_ratio_name)
        os.unlink(self.fasta_file_for_consensus_unequal_length_name)
        
        # delete = False in creation of tempfiles
        # as they need to be opened in functions
        # hence the need to remove them

    def test_select_unique_rand_bcs(self):
        fasta_seqs_of_rand_bcs = self.fasta_seqs_of_rand_bcs
        actual = select_unique_rand_bcs(fasta_seqs_of_rand_bcs, 0.86)
        expected = set(['ATTGCATTGCATTGCATTGC', 'ATTGCTTATTGCATTGCTTT'])
        self.assertEqual(actual, expected)

    def test_get_consensus(self):
        fasta_seqs_for_consensus = self.fasta_seqs_for_consensus
        temp_file = self.fasta_file_for_consensus
        temp_file.write(fasta_seqs_for_consensus_tie_G_C)
        # at the last position, G and C have the same frequency
        # therefore the function is expected to return
        # consensus sequence with G, which is present in seq 
        # that appears max times. (10, 10) while C appreared
        # in sequence that have count: (9, 6, 5)
        # If there is still a tie, the function will return
        # the base that appeared first.
        # This method is just for a consistent way
        # to resolve ties        
        temp_file_name = temp_file.name
        temp_file.close()
        temp_file = open(temp_file_name, 'r')
        actual = get_consensus(temp_file, 2)
        expected = 'ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG'
        self.assertEqual(actual, expected)

        temp_file = self.fasta_file_for_consensus_unequal_length
        temp_file.write(fasta_seqs_for_consensus_unequal_length)
        # Expected to throw SeqLengthMismatchError
        temp_file_name = temp_file.name
        temp_file.close()
        temp_file = open(temp_file_name, 'r')
        self.assertRaises(SeqLengthMismatchError)

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
        
    def test_extract_primers(self):
        fasta_seq_for_primer = self.fasta_seq_for_primer
        possible_primers = self.possible_primers
        actual = extract_primer(fasta_seq_for_primer, possible_primers)
        expected = ('A', 'ATGC', 'CCCC')
        self.assertEqual(actual, expected)        

fasta_seqs_for_uclust = """>1abc|1
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
>2abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>3abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>4abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>5abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>6abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>7abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>8abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>9abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
"""

fasta_seqs_for_consensus = """>1id1|1
ATGCATGG
>2id2|14
ATGCATGC
"""

fasta_seqs_for_consensus_unequal_length = """>1id1|1
ATGCATGG
>2id2|14
ATGCATGCT
"""

fasta_seqs_for_consensus_tie_G_C = """>abc1|10
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
>abc1|9
ATTTTATGGGCGGCGCGCCGCGCGCGCATTATATATATATAGCGCGCGCGCGCGC
>abc1|5
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGC
>abc1|10
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
>abc1|6
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGC
"""

fasta_seqs_of_rand_bcs = ['ATTGCATTGCATTGCATTGC', 'ATTGCATTGCATTGCATTGC', 'ATTGCATTGCATTGCATTG', 'ATTGCTTATTGCATTGCTTT']

# run tests if called from command line
if __name__ == '__main__':
    main()
