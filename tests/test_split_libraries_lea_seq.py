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
from unittest import TestCase, main
from qiime.util import get_qiime_temp_dir
from qiime.split_libraries_lea_seq import (get_cluster_ratio, get_consensus,
                                           get_LEA_seq_consensus_seqs,
                                           select_unique_rand_bcs,
                                           extract_primer,
                                           SeqLengthMismatchError)
from skbio.util.misc import remove_files
import os


class WorkflowTests(TestCase):

    def setUp(self):
        """setup the test values"""
        # define test data
        self.fasta_seqs_of_rand_bcs = fasta_seqs_of_rand_bcs
        self.fasta_seqs_for_uclust = fasta_seqs_for_uclust
        self.fasta_seqs_for_consensus = fasta_seqs_for_consensus
        self.fwd_read_data = fwd_read_data
        self.rev_read_data = rev_read_data
        self.mapping_data = mapping_data
        self.fasta_seq_for_primer = fasta_seq_for_primer
        self.possible_primers = possible_primers
        self.fasta_seqs_for_consensus_tie_G_C = \
            fasta_seqs_for_consensus_tie_G_C
        self.fasta_seqs_for_consensus_unequal_length = \
            fasta_seqs_for_consensus_unequal_length
        self.min_difference_in_clusters = min_difference_in_clusters

        # open temporary files needed
        temp_dir = get_qiime_temp_dir()
        self.fasta_file_for_consensus_tie_G_C = tempfile.NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=temp_dir)
        self.fasta_file_for_cluster_ratio = tempfile.NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=temp_dir)
        self.fasta_file_for_consensus_unequal_length = \
            tempfile.NamedTemporaryFile(delete=False,
                                        mode='w',
                                        dir=temp_dir)
        self.mapping_fp = tempfile.NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=temp_dir)
        self.rev_read_fp = tempfile.NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=temp_dir)
        self.fwd_read_fp = tempfile.NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=temp_dir)
        self.log_file = tempfile.NamedTemporaryFile(mode='w', dir=temp_dir)

        # define temp dir
        self.temp_dir = temp_dir

        # write test data in temp files
        self.mapping_fp.write(self.mapping_data)
        self.fwd_read_fp.write(self.fwd_read_data)
        self.rev_read_fp.write(self.rev_read_data)
        self.fasta_file_for_cluster_ratio.write(self.fasta_seqs_for_uclust)
        self.fasta_file_for_consensus_tie_G_C.write(
            self.fasta_seqs_for_consensus_tie_G_C)
        self.fasta_file_for_consensus_unequal_length.write(
            self.fasta_seqs_for_consensus_unequal_length)
        # self.log_file does not need any data

        # save names of the files
        # so it can be used for system calls
        # after they are closed
        self.fasta_file_for_consensus_tie_G_C_name = \
            self.fasta_file_for_consensus_tie_G_C.name
        self.fasta_file_for_consensus_unequal_length_name = \
            self.fasta_file_for_consensus_unequal_length.name
        self.fasta_file_for_cluster_ratio_name = \
            self.fasta_file_for_cluster_ratio.name
        self.mapping_fp_name = self.mapping_fp.name
        self.rev_read_fp_name = self.rev_read_fp.name
        self.fwd_read_fp_name = self.fwd_read_fp.name

        self.sequence_read_fps = (self.fwd_read_fp.name, self.rev_read_fp.name)
        # close and open the files that are needed to be in read mode
        self.fasta_file_for_consensus_tie_G_C.close()
        self.fasta_file_for_consensus_unequal_length.close()
        self.fasta_file_for_consensus_tie_G_C = open(
            self.fasta_file_for_consensus_tie_G_C_name,
            'r')
        self.fasta_file_for_consensus_unequal_length = open(
            self.fasta_file_for_consensus_unequal_length_name,
            'r')

        # close files that don't need to be open
        # 1. (some of which are only needed for system calls)
        # 2. test_get_LEA_seq_consensus_seqs is the main library function
        # which only requires names, and it will open the files
        self.fasta_file_for_cluster_ratio.close()
        self.mapping_fp.close()
        self.rev_read_fp.close()
        self.fwd_read_fp.close()

    def tearDown(self):
        """remove all the files after completing tests """
        # delete = False in creation of tempfiles
        # as they need to be opened in functions
        # hence the need to remove them
        files_to_be_removed = list()
        files_to_be_removed.append(self.fasta_file_for_consensus_tie_G_C_name)
        files_to_be_removed.append(
            self.fasta_file_for_consensus_unequal_length_name)
        files_to_be_removed.append(self.fasta_file_for_cluster_ratio_name)
        files_to_be_removed.append(self.mapping_fp_name)
        files_to_be_removed.append(self.rev_read_fp_name)
        files_to_be_removed.append(self.fwd_read_fp_name)
        remove_files(files_to_be_removed)
        self.log_file.close()  # will be deleted

    def test_select_unique_rand_bcs(self):
        fasta_seqs_of_rand_bcs = self.fasta_seqs_of_rand_bcs
        actual = select_unique_rand_bcs(fasta_seqs_of_rand_bcs, 0.86)
        expected = set(['ATTGCATTGCATTGCATTGC', 'ATTGCTTATTGCATTGCTTT'])
        self.assertEqual(actual, expected)

    def test_get_consensus(self):
        fasta_file_for_consensus_unequal_length = \
            self.fasta_file_for_consensus_unequal_length
        fasta_file_for_consensus_tie_G_C = \
            self.fasta_file_for_consensus_tie_G_C

        actual = get_consensus(fasta_file_for_consensus_tie_G_C, 2)
        # at the last position, G and C have the same frequency
        # therefore the function is expected to return
        # consensus sequence with G, which is present in seq
        # that appears max times. (10, 10) while C appreared
        # in sequence that have count: (9, 6, 5)
        # If there is still a tie, the function will return
        # the base that appeared first.
        # This method is just for a consistent way
        # to resolve ties
        expected = 'ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG'
        self.assertEqual(actual, expected)

        # Sequences having unequal length:
        self.assertRaises(
            SeqLengthMismatchError,
            get_consensus,
            fasta_file_for_consensus_unequal_length,
            2)

    def test_get_cluster_ratio(self):
        min_difference_in_clusters = self.min_difference_in_clusters
        fasta_file_for_cluster_ratio_name = \
            self.fasta_file_for_cluster_ratio_name
        actual = get_cluster_ratio(
            fasta_file_for_cluster_ratio_name,
            min_difference_in_clusters)
        expected = 0.125
        self.assertEqual(actual, expected)

    def test_extract_primers(self):
        fasta_seq_for_primer = self.fasta_seq_for_primer
        possible_primers = self.possible_primers
        actual = extract_primer(fasta_seq_for_primer, possible_primers)
        expected = ('A', 'ATGC', 'CCCC')
        self.assertEqual(actual, expected)

    def test_get_LEA_seq_consensus_seqs(self):
        sequence_read_fps = self.sequence_read_fps
        mapping_fp_name = self.mapping_fp_name
        temp_dir = self.temp_dir
        barcode_type = int(7)
        barcode_len = 7
        barcode_correction_fn = None
        max_barcode_errors = 1.5
        min_consensus = 0.66
        max_cluster_ratio = 2.5
        min_difference_in_bcs = 0.86
        log_file = self.log_file
        fwd_length = 19
        rev_length = 19
        min_reads_per_random_bc = 1
        min_difference_in_clusters = self.min_difference_in_clusters
        function_call = get_LEA_seq_consensus_seqs(sequence_read_fps,
                                                   mapping_fp_name, temp_dir,
                                                   barcode_type, barcode_len,
                                                   barcode_correction_fn,
                                                   max_barcode_errors,
                                                   min_consensus,
                                                   max_cluster_ratio,
                                                   min_difference_in_bcs,
                                                   log_file, fwd_length,
                                                   rev_length,
                                                   min_reads_per_random_bc,
                                                   min_difference_in_clusters)

        actual = function_call['Sample1']['AGCTACGAGCTATTGC']
        expected = 'AAAAAAAAAAAAAAAAAAA^AAAAAAAAAAAAAAAAAA'
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

fasta_seqs_of_rand_bcs = [
    'ATTGCATTGCATTGCATTGC',
    'ATTGCATTGCATTGCATTGC',
    'ATTGCATTGCATTGCATTG',
    'ATTGCTTATTGCATTGCTTT']


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
Sample1	CCGGCAG	AGAGTTTGATCCTGGCTCAG	GGGCCGTGTCTCAGT	Sample1	description
"""

min_difference_in_clusters = 0.98
fasta_seq_for_primer = 'AATGCCCCC'
possible_primers = ['ATGC', 'ATTT']

# run tests if called from command line
if __name__ == '__main__':
    main()
