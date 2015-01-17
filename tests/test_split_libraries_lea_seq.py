#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"

from tempfile import NamedTemporaryFile
from unittest import TestCase, main

from skbio.util import remove_files
from qiime.util import get_qiime_temp_dir
from qiime.split_libraries_lea_seq import (get_cluster_ratio, get_consensus,
                                           get_LEA_seq_consensus_seqs,
                                           select_unique_rand_bcs,
                                           extract_primer,
                                           format_lea_seq_log,
                                           process_mapping_file,
                                           check_barcodes,
                                           get_consensus_seqs_lookup,
                                           read_fwd_rev_read,
                                           InvalidGolayBarcodeError,
                                           BarcodeLenMismatchError,
                                           SeqLengthMismatchError,
                                           LowConsensusScoreError,
                                           PrimerMismatchError)


class WorkflowTests(TestCase):

    def setUp(self):
        """setup the test values"""
        # define test data
        self.fasta_seqs_of_rand_bcs = fasta_seqs_of_rand_bcs
        self.fasta_seqs_for_cluster_ratio = fasta_seqs_for_cluster_ratio
        self.fasta_seqs_for_consensus = fasta_seqs_for_consensus
        self.fwd_read_data = fwd_read_data.split()
        self.rev_read_data = rev_read_data.split()
        self.mapping_data = mapping_data
        self.fasta_seq_for_primer = fasta_seq_for_primer
        self.possible_primers = possible_primers
        self.fasta_seqs_for_consensus_tie_G_C = \
            fasta_seqs_for_consensus_tie_G_C
        self.fasta_seqs_for_consensus_unequal_length = \
            fasta_seqs_for_consensus_unequal_length
        self.min_difference_in_clusters = min_difference_in_clusters
        self.temp_dir = get_qiime_temp_dir()
        self.mapping_fp = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.mapping_fp.write(self.mapping_data)
        self.mapping_fp_name = self.mapping_fp.name
        self.mapping_fp.close()
        self.mapping_fp = open(self.mapping_fp_name, 'r')
        self.seqs_with_no_consensus = seqs_with_no_consensus
        self.false_primers = false_primers
        self.barcode_len = barcode_len
        self.barcode_correction_fn = barcode_correction_fn
        self.max_barcode_errors = max_barcode_errors
        self.fwd_length = fwd_length
        self.rev_length = fwd_length
        self.bc_to_sid = bc_to_sid
        self.bc_to_fwd_primers = bc_to_fwd_primers
        self.bc_to_rev_primers = bc_to_rev_primers
        self.min_difference_in_bcs = min_difference_in_bcs
        self.min_reads_per_random_bc = min_reads_per_random_bc
        self.max_cluster_ratio = max_cluster_ratio

    def tearDown(self):
        """remove all the files after completing tests """
        self.mapping_fp.close()
        remove_files([self.mapping_fp_name])

    def test_select_unique_rand_bcs(self):
        actual = select_unique_rand_bcs(self.fasta_seqs_of_rand_bcs, 0.86)
        expected = set(['ATTGCATTGCATTGCATTGC', 'ATTGCTTATTGCATTGCTTT'])
        self.assertEqual(actual, expected)

    def test_get_consensus(self):
        actual = get_consensus(self.fasta_seqs_for_consensus_tie_G_C, 2)
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
        with self.assertRaises(SeqLengthMismatchError):
            get_consensus(self.fasta_seqs_for_consensus_unequal_length, 2)

        seqs_with_no_consensus = self.seqs_with_no_consensus
        with self.assertRaises(LowConsensusScoreError):
            get_consensus(seqs_with_no_consensus, 6.6)

    def test_get_cluster_ratio(self):
        actual = get_cluster_ratio(
            self.fasta_seqs_for_cluster_ratio,
            self.min_difference_in_clusters)
        expected = 0.125
        self.assertEqual(actual, expected)

    def test_extract_primers(self):
        actual = extract_primer(
            self.fasta_seq_for_primer, self.possible_primers)
        expected = ('A', 'ATGC', 'CCCC')
        self.assertEqual(actual, expected)
        with self.assertRaises(PrimerMismatchError):
            extract_primer(self.fasta_seq_for_primer, self.false_primers)

    def test_get_LEA_seq_consensus_seqs(self):
        barcode_type = int(7)
        barcode_len = 7
        barcode_correction_fn = None
        max_barcode_errors = 1.5
        min_consensus = 0.66
        max_cluster_ratio = 2.5
        min_difference_in_bcs = 0.86
        fwd_length = 19
        rev_length = 19
        min_reads_per_random_bc = 1
        min_diff_in_clusters = self.min_difference_in_clusters
        barcode_column = 'BarcodeSequence'
        reverse_primer_column = 'ReversePrimer'

        function_call, _ = get_LEA_seq_consensus_seqs(self.fwd_read_data,
                                                      self.rev_read_data,
                                                      self.mapping_fp,
                                                      self.temp_dir,
                                                      barcode_type,
                                                      barcode_len,
                                                      barcode_correction_fn,
                                                      max_barcode_errors,
                                                      min_consensus,
                                                      max_cluster_ratio,
                                                      min_difference_in_bcs,
                                                      fwd_length,
                                                      rev_length,
                                                      min_reads_per_random_bc,
                                                      min_diff_in_clusters,
                                                      barcode_column,
                                                      reverse_primer_column)

        actual = function_call['Sample1']['AGCTACGAGCTATTGC']
        expected = 'AAAAAAAAAAAAAAAAAAA^AAAAAAAAAAAAAAAAAA'
        self.assertEqual(actual, expected)

    def test_format_lea_seq_log(self):
        actual = format_lea_seq_log(1, 2, 3, 4, 5, 6)
        expected = """Quality filter results
Total number of input sequences: 1
Barcode not in mapping file: 3
Sequence shorter than threshold: 5
Barcode errors exceeds limit: 2
Primer mismatch count: 4

Total number seqs written: 6"""
        self.assertEqual(actual, expected)

    def test_process_mapping_file(self):
        barcode_type = int(7)
        barcode_len = 7
        barcode_column = 'BarcodeSequence'
        reverse_primer_column = 'ReversePrimer'

        actual = process_mapping_file(self.mapping_fp,
                                      barcode_len,
                                      barcode_type,
                                      barcode_column,
                                      reverse_primer_column)
        bc_to_sid = ({'CCGGCAG': 'Sample1'},)
        bc_to_fwd_primers = ({'CCGGCAG': {'AGAGTTTGATCCTGGCTCAG': 20}},)
        bc_to_rev_primers = ({'CCGGCAG': ['GGGCCGTGTCTCAGT']},)
        expected = bc_to_sid + bc_to_fwd_primers + bc_to_rev_primers
        self.assertEqual(actual, expected)

    def test_check_barcodes(self):
        barcode_type = 'golay_12'
        barcode_len = 7
        bc_to_sid = {'CCGGCAG': 'Sample1'}
        with self.assertRaises(InvalidGolayBarcodeError):
            check_barcodes(bc_to_sid, barcode_len, barcode_type)
        barcode_len = 1
        with self.assertRaises(BarcodeLenMismatchError):
            check_barcodes(bc_to_sid, barcode_len, barcode_type)

    def test_read_fwd_rev_read(self):
        expected_seqs_kept = 1
        function_call = read_fwd_rev_read(self.fwd_read_data,
                                          self.rev_read_data,
                                          self.bc_to_sid,
                                          self.barcode_len,
                                          self.barcode_correction_fn,
                                          self.bc_to_fwd_primers,
                                          self.bc_to_rev_primers,
                                          self.max_barcode_errors,
                                          self.fwd_length,
                                          self.rev_length)
        actual_seqs_kept = function_call[-1]
        self.assertEqual(actual_seqs_kept, expected_seqs_kept)

    def test_get_consensus_seqs_lookup(self):

        fn_call_fwd_rev_read = read_fwd_rev_read(self.fwd_read_data,
                                                 self.rev_read_data,
                                                 self.bc_to_sid,
                                                 self.barcode_len,
                                                 self.barcode_correction_fn,
                                                 self.bc_to_fwd_primers,
                                                 self.bc_to_rev_primers,
                                                 self.max_barcode_errors,
                                                 self.fwd_length,
                                                 self.rev_length)

        random_bc_lookup = fn_call_fwd_rev_read[0]
        random_bc_reads = fn_call_fwd_rev_read[1]
        random_bcs = fn_call_fwd_rev_read[2]
        min_difference_bcs = self.min_difference_in_bcs
        min_diff_clusters = self.min_difference_in_clusters
        min_reads_rand_bc = self.min_reads_per_random_bc
        max_cluster_ratio = self.max_cluster_ratio
        output_dir = self.temp_dir

        fn_call_get_consensus = get_consensus_seqs_lookup(random_bc_lookup,
                                                          random_bc_reads,
                                                          random_bcs,
                                                          min_difference_bcs,
                                                          min_reads_rand_bc,
                                                          output_dir,
                                                          min_diff_clusters,
                                                          max_cluster_ratio,
                                                          min_consensus)

        actual = fn_call_get_consensus['Sample1']['AGCTACGAGCTATTGC']
        expected = 'AAAAAAAAAAAAAAAAAAA^AAAAAAAAAAAAAAAAAA'
        self.assertEqual(actual, expected)


fasta_seqs_for_cluster_ratio = """>1abc|1
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
mapping_data = """"#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	Description
Sample1	CCGGCAG	AGAGTTTGATCCTGGCTCAG	GGGCCGTGTCTCAGT	Sample1	description"""

barcode_type = int(7)
barcode_len = 7
barcode_correction_fn = None
max_barcode_errors = 1.5
min_consensus = 0.66
max_cluster_ratio = 2.5
min_difference_in_bcs = 0.86
fwd_length = 19
rev_length = 19
min_reads_per_random_bc = 1
barcode_column = 'BarcodeSequence'
reverse_primer_column = 'ReversePrimer'
bc_to_sid = {'CCGGCAG': 'Sample1'}
bc_to_fwd_primers = {'CCGGCAG': {'AGAGTTTGATCCTGGCTCAG': 20}}
bc_to_rev_primers = {'CCGGCAG': ['GGGCCGTGTCTCAGT']}


seqs_with_no_consensus = """>id1|1
ATGC
>id2|1
TGCA
>id3|1
GCAT
>id4|1
CATG"""


min_difference_in_clusters = 0.98
fasta_seq_for_primer = 'AATGCCCCC'
possible_primers = ['ATGC', 'ATTT']
false_primers = ['AAAA']
# run tests if called from command line
if __name__ == '__main__':
    main()
