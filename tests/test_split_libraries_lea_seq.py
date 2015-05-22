#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
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
        self.temp_dir = get_qiime_temp_dir()
        self.fasta_seqs_of_rand_bcs = fasta_seqs_of_rand_bcs
        self.fasta_seqs_for_cluster_ratio = fasta_seqs_for_cluster_ratio
        self.fasta_seqs_for_consensus = fasta_seqs_for_consensus
        self.fwd_read_data = fwd_read_data
        self.rev_read_data = rev_read_data
        self.get_cons_fwd_read_data = get_cons_fwd_read_data
        self.get_cons_rev_read_data = get_cons_rev_read_data

        self.fwd_read_fh = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.fwd_read_fh_name = self.fwd_read_fh.name
        self.fwd_read_fh.write(self.fwd_read_data)
        self.fwd_read_fh.close()
        self.fwd_read_fh = open(self.fwd_read_fh_name, 'r')
        self.rev_read_fh = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.rev_read_fh_name = self.rev_read_fh.name
        self.rev_read_fh.write(self.rev_read_data)
        self.rev_read_fh.close()
        self.rev_read_fh = open(self.rev_read_fh_name, 'r')

        self.get_cons_fwd_read_fh = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.get_cons_fwd_read_fh_name = self.get_cons_fwd_read_fh.name
        self.get_cons_fwd_read_fh.write(self.get_cons_fwd_read_data)
        self.get_cons_fwd_read_fh.close()
        self.get_cons_fwd_read_fh = open(self.get_cons_fwd_read_fh_name, 'r')
        self.get_cons_rev_read_fh = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.get_cons_rev_read_fh_name = self.get_cons_rev_read_fh.name
        self.get_cons_rev_read_fh.write(self.get_cons_rev_read_data)
        self.get_cons_rev_read_fh.close()
        self.get_cons_rev_read_fh = open(self.get_cons_rev_read_fh_name, 'r')

        self.mapping_data = mapping_data
        self.get_cons_mapping_data = get_cons_mapping_data
        self.fasta_seq_for_primer = fasta_seq_for_primer
        self.possible_primers = possible_primers

        self.fasta_seqs_for_consensus_tie_G_C = \
            fasta_seqs_for_consensus_tie_G_C
        self.fasta_seqs_for_consensus_unequal_length = \
            fasta_seqs_for_consensus_unequal_length
        self.seqs_with_no_consensus = seqs_with_no_consensus

        self.fasta_file_for_consensus_tie_G_C = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.fasta_file_for_consensus_tie_G_C_name = \
            self.fasta_file_for_consensus_tie_G_C.name
        self.fasta_file_for_consensus_tie_G_C.write(
            self.fasta_seqs_for_consensus_tie_G_C)
        self.fasta_file_for_consensus_tie_G_C.close()
        self.fasta_file_for_consensus_tie_G_C = open(
            self.fasta_file_for_consensus_tie_G_C_name, 'r')

        self.fasta_file_for_consensus_unequal_length = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.fasta_file_for_consensus_unequal_length_name = \
            self.fasta_file_for_consensus_unequal_length.name
        self.fasta_file_for_consensus_unequal_length.write(
            self.fasta_seqs_for_consensus_unequal_length)
        self.fasta_file_for_consensus_unequal_length.close()
        self.fasta_file_for_consensus_unequal_length = open(
            self.fasta_file_for_consensus_unequal_length_name, 'r')

        self.fasta_file_no_consensus = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.fasta_file_no_consensus_name = self.fasta_file_no_consensus.name
        self.fasta_file_no_consensus.write(self.seqs_with_no_consensus)
        self.fasta_file_no_consensus.close()
        self.fasta_file_no_consensus = open(
            self.fasta_file_no_consensus_name, 'r')

        self.min_difference_in_clusters = min_difference_in_clusters

        self.mapping_fp = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.mapping_fp.write(self.mapping_data)
        self.mapping_fp_name = self.mapping_fp.name
        self.mapping_fp.close()
        self.mapping_fp = open(self.mapping_fp_name, 'r')

        self.get_cons_mapping_fp = NamedTemporaryFile(
            delete=False,
            mode='w',
            dir=self.temp_dir)
        self.get_cons_mapping_fp.write(self.get_cons_mapping_data)
        self.get_cons_mapping_fp_name = self.get_cons_mapping_fp.name
        self.get_cons_mapping_fp.close()
        self.get_cons_mapping_fp = open(self.get_cons_mapping_fp_name, 'r')

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
        self.fasta_file_no_consensus.close()
        self.fasta_file_for_consensus_tie_G_C.close()
        self.fasta_file_for_consensus_unequal_length.close()
        remove_files([self.mapping_fp_name,
                      self.fasta_file_no_consensus_name,
                      self.fasta_file_for_consensus_tie_G_C_name,
                      self.fasta_file_for_consensus_unequal_length_name,
                      self.fwd_read_fh_name, self.rev_read_fh_name])

    def test_select_unique_rand_bcs(self):
        actual = select_unique_rand_bcs(self.fasta_seqs_of_rand_bcs, 0.86)
        expected = set(['ATTGCATTGCATTGCATTGC', 'ATTGCTTATTGCATTGCTTT'])
        self.assertEqual(actual, expected)

    def test_get_consensus(self):
        actual = get_consensus(self.fasta_file_for_consensus_tie_G_C, 2)
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
            get_consensus(self.fasta_file_for_consensus_unequal_length, 2)

        fasta_file_no_consensus = self.fasta_file_no_consensus
        with self.assertRaises(LowConsensusScoreError):
            get_consensus(fasta_file_no_consensus, 6.6)

    def test_get_cluster_ratio(self):
        actual = get_cluster_ratio(
            self.fasta_seqs_for_cluster_ratio,
            self.min_difference_in_clusters)
        expected = 2.5
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

        function_call, _ = get_LEA_seq_consensus_seqs(self.fwd_read_fh,
                                                      self.rev_read_fh,
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
        # this call tests the second condition of if loop
        # in the function get_consensus_seq_lookup
        # i.e. select the majority sequence, as the cluster ratio
        # between max_cluster/second_best_cluster in the fwd_read_data
        # (and rev_read_data) is 3/1 > 2.5,
        # so the function get_consensus will not be called

        fn_call, _ = get_LEA_seq_consensus_seqs(self.get_cons_fwd_read_fh,
                                                self.get_cons_rev_read_fh,
                                                self.get_cons_mapping_fp,
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

        get_cons_actual = fn_call['Sample1']['AGCTACGAGCTATTGC']
        get_cons_expected = 'AAAAAAAAAACAAAAAAAA^AAAAAAAAAATAAAAATA'
        self.assertEqual(get_cons_actual, get_cons_expected)
        # this call tests the third condition of if loop
        # in the function get_consensus_seq_lookup.
        # i.e. calls the get_consensus function, as the cluster ratio
        # between max_cluster/second_best_cluster in the get_cons_fwd_read_data
        # (and get_cons_rev_read_data) is 2/1 ( < 2.5)
        # so the majority sequence will not be selected

        get_cons_actual = fn_call['Sample2']['AGCTACGCATCAAGGG']
        get_cons_expected = 'AAAAAAAAAATAAAAAAAA^TTAAAAAAAAAAAAGAAAA'
        self.assertEqual(get_cons_actual, get_cons_expected)

        self.assertFalse(len(fn_call) <= 1,
                         msg="The get_consensus_seqs_lookup function "
                         "has returned early, without completing "
                         "the three 'for' loops.")

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
        expected_seqs_kept = 4
        function_call = read_fwd_rev_read(self.fwd_read_fh,
                                          self.rev_read_fh,
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

        fn_call_fwd_rev_read = read_fwd_rev_read(self.fwd_read_fh,
                                                 self.rev_read_fh,
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
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>4abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>5abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>6abc|1
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>7abc|1
ATTTTATTTTATTTTTATTTATTATATATTATATATATATAGCGCGCGCGCGCGG
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>8abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>9abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>10abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>11abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>12abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>13abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
>14abc|1
GGTCGGTCGTGCGTGCTCGTCGTGCTCGTCGTCGTCGCTCGTCGTCGCTGCTCTC
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@2____
AGCTACGAGCTATTGCAGAGTTTGATCCTGGCTCAGAAAAAAAAAAAAAAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@3____
AGCTACGAGCTATTGCAGAGTTTGATCCTGGCTCAGAAAAAAAAAAAAAAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@4____
AGCTACGAGCTATTGCAGAGTTTGATCCTGGCTCAGAAAAAAAAAAATTAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
rev_read_data = """@1____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAAAAAAAAAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@2____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAAAAAAAAAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@3____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAAAAAAAAAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@4____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAAAAAAACCA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
mapping_data = """"#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	Description
Sample1	CCGGCAG	AGAGTTTGATCCTGGCTCAG	GGGCCGTGTCTCAGT	Sample1	description"""

get_cons_fwd_read_data = """@1____
AGCTACGCATCAAGGGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAGAAAAAAAACCAACAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@2____
AGCTACGCATCAAGGGTTTTTTTTTTTTTTTTTTTTAAAAAAAAAATAAAAAAAACCAACAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@3____
AGCTACGAGCTATTGCAGAGTTTGATCCTGGCTCAGAAAAAAAAAACAAAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@4____
AGCTACGAGCTATTGCAGAGTTTGATCCTGGCTCAGAAAAAAAAAACAAAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@5____
AGCTACGAGCTATTGCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAGAAAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@6____
AGCTACGAGCTATTGCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAATAAAAAAAACCGGCAG
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""

get_cons_rev_read_data = """@1____
CCAACAGAGCTACGAGCTATTTTTTTTTTTTTTTTTAAAAAAAAAAAAGAAAAAAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@2____
CCAACAGAGCTACGAGCTATTTTTTTTTTTTTTTTTAAAAAAAAAAAAGAAAAAAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@3____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAATAAAAACA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@4____
CCGGCAGAGCTACGAGCTATTGCGGGCCGTGTCTCAGTAAAAAAAAAAAAAAAATA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@5____
CCGGCAGAGCTACGAGCTATTTTTTTTTTTTTTTTAAAAAAAAAAAATAAAAACAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
@6____
CCGGCAGAGCTACGAGCTATTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAATAA
+
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
"""
get_cons_mapping_data = """"#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	Description
Sample1	CCGGCAG	AGAGTTTGATCCTGGCTCAG	GGGCCGTGTCTCAGT	Sample1	description
Sample2	CCAACAG	TTTTTTTTTTTTTTTTTTTT	TTTTTTTTTTTTTTT	Sample2	description"""

# breakdown of get_cons_fwd_read_data = """@1____
# for testing:
# AGCTACGCATCAAGGG random barcode sequence 1-16
# AGAGTTTGATCCTGGCTCAG Linker Primer sequence 17 - 36
# AAAAAAAAAAGAAAAAAAA sequence 37 - 55
# CCGGCAG BarcodeSequence 56 - 63

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
