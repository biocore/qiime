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

        self.fwd_read_data = fwd_read_data
        self.rev_read_data = rev_read_data
        self.mapping_data = mapping_data
        self.mapping_fp = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=temp_dir)
        self.rev_read_fp = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=temp_dir)
        self.fwd_read_fp = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=temp_dir)
        self.mapping_fp.write(self.mapping_data)
        self.fwd_read_fp.write(self.fwd_read_data)
        self.rev_read_fp.write(self.rev_read_data)
        self.sequence_read_fps = (self.fwd_read_fp.name ,self.rev_read_fp.name)
        self.log_file = tempfile.NamedTemporaryFile(mode='w', dir=temp_dir)
        self.mapping_fp_name = self.mapping_fp.name
        self.rev_read_fp_name = self.rev_read_fp.name
        self.fwd_read_fp_name = self.fwd_read_fp.name
        self.mapping_fp.close()
        self.rev_read_fp.close()
        self.fwd_read_fp.close()
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
        self.log_file.close()
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

    def test_get_LEA_seq_consensus_seqs(self):
        fwd_read_data = self.fwd_read_data
        rev_read_data = self.rev_read_data
        mapping_data = self.mapping_data
        temp_dir = get_qiime_temp_dir()
        temp_dir = temp_dir 
        barcode_type = int(7)
        barcode_correction_fn = None
        max_barcode_errors = 1.5
        min_consensus = 0.66
        max_cluster_ratio = 2.5
        min_difference_in_bcs = 0.86
        log_file = self.log_file
        function_call = get_LEA_seq_consensus_seqs(self.sequence_read_fps, self.mapping_fp.name,
                                           temp_dir, barcode_type, barcode_correction_fn,
                                           max_barcode_errors, min_consensus,
                                           max_cluster_ratio, min_difference_in_bcs, log_file)
        actual = function_call['Sample1']['AGCTACGAGCTATTGC']
        expected = 'AAAAAAAAAAAAAAAAAAA, AAAAAAAAAAAAAAAAAA'
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


fwd_read_data = """@1____
AGCTACGAGCTATTGCAGAGTTTGATCCTGGCTCAGAAAAAAAAAAAAAAAAAAACCGGCAG
+
$
"""
rev_read_data = """@1____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAAAAAAAAAA
+
$
"""
mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	Description
Sample1	CCGGCAG	AGAGTTTGATCCTGGCTCAG	GGGCCGTGTCTCAGT	Sample1 description
"""


# run tests if called from command line
if __name__ == '__main__':
    main()
