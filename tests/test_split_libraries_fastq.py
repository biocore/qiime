#!/usr/bin/env python
# File created on 05 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import numpy as np

from unittest import TestCase, main
from qiime.split_libraries_fastq import (
    process_fastq_single_end_read_file,
    quality_filter_sequence,
    bad_chars_from_threshold,
    get_illumina_qual_chars,
    quality_filter_sequence,
    FastqParseError,
    check_header_match_pre180,
    check_header_match_180_or_later,
    correct_barcode,
    process_fastq_single_end_read_file_no_barcode
)
from qiime.golay import decode_golay_12

# yes, this is terrible but so is manually converting all the qual scores.
from skbio.parse.sequences.fastq import _ascii_to_phred64, _ascii_to_phred33

class FakeFile(object):

    def __init__(self, d=""):
        self.s = d

    def write(self, s):
        self.s += s

    def close(self):
        pass


class SplitLibrariesFastqTests(TestCase):

    """ """

    def setUp(self):

        self.fastq1 = fastq1.split('\n')
        self.barcode_fastq1 = barcode_fastq1.split('\n')
        self.fastq2 = fastq2.split('\n')
        self.barcode_fastq2 = barcode_fastq2.split('\n')
        self.fastq1_expected_no_qual_unassigned = fastq1_expected_no_qual_unassigned
        self.fastq1_expected_default = fastq1_expected_default
        self.fastq2_expected_default = fastq2_expected_default
        self.fastq1_expected_single_barcode = fastq1_expected_single_barcode
        self.barcode_map1 = barcode_map1

    def test_correct_barcode_exact_match(self):
        """correct_barcode functions as expected w exact match"""
        barcode = "GGAGACAAGGGA"
        barcode_to_sample_id = {
            "GGAGACAAGGGA": "s1",
            "ACACCTGGTGAT": "s2"}
        correction_fn = None

        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (0, barcode, False, 's1')
        self.assertEqual(actual, expected)

        correction_fn = decode_golay_12
        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (0, barcode, False, 's1')
        self.assertEqual(actual, expected)

    def test_correct_barcode_no_error_correction(self):
        """correct_barcode functions as expected w no error correction"""
        barcode = "GGAGACAAGGGT"
        barcode_to_sample_id = {
            "GGAGACAAGGGA": "s1",
            "ACACCTGGTGAT": "s2"}
        correction_fn = None

        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (0, barcode, False, None)
        self.assertEqual(actual, expected)

        # barcode contains N
        barcode = "CCAGTGTANGCA"
        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (0, "CCAGTGTANGCA", False, None)
        self.assertEqual(actual, expected)

    def test_correct_barcode_golay_correction(self):
        """correct_barcode functions as expected w golay correction"""
        barcode = "GGAGACAAGGGT"
        barcode_to_sample_id = {
            "GGAGACAAGGGA": "s1",
            "ACACCTGGTGAT": "s2"}
        correction_fn = decode_golay_12

        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (1, "GGAGACAAGGGA", True, "s1")
        self.assertEqual(actual, expected)

        barcode = "ACACCTGGTGAC"
        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (1, "ACACCTGGTGAT", True, "s2")
        self.assertEqual(actual, expected)

        # valid code, but not in barcode_to_sample_id map
        barcode = "CCAGTGTATGCA"
        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (0, "CCAGTGTATGCA", True, None)
        self.assertEqual(actual, expected)

        # invalid code, corrected not in barcode_to_sample_id map
        barcode = "CCTGTGTATGCA"
        actual = correct_barcode(barcode, barcode_to_sample_id, correction_fn)
        expected = (1, "CCAGTGTATGCA", True, None)
        self.assertEqual(actual, expected)

    def test_process_fastq_single_end_read_file(self):
        """process_fastq_single_end_read_file functions as expected w no qual filter
        """
        actual = process_fastq_single_end_read_file(self.fastq1,
                                                    self.barcode_fastq1,
                                                    self.barcode_map1,
                                                    store_unassigned=True,
                                                    max_bad_run_length=1000,
                                                    phred_quality_threshold=None,
                                                    min_per_read_length_fraction=0.,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=1000,
                                                    start_seq_id=0,
                                                    filter_bad_illumina_qual_digit=False)
        actual = list(actual)
        expected = self.fastq1_expected_no_qual_unassigned
        self.assertEqual(len(actual), len(expected))
        for i in range(len(expected)):
            np.testing.assert_equal(actual[i], expected[i])

    def test_process_fastq_single_end_read_file_w_defaults(self):
        """process_fastq_single_end_read_file functions as expected w default filters
        """
        actual = process_fastq_single_end_read_file(self.fastq1,
                                                    self.barcode_fastq1,
                                                    self.barcode_map1,
                                                    min_per_read_length_fraction=0.45)
        actual = list(actual)
        expected = self.fastq1_expected_default
        self.assertEqual(len(actual), len(expected))
        for i in range(len(expected)):
            np.testing.assert_equal(actual[i], expected[i])

    def test_process_fastq_single_end_read_file_no_barcode(self):
        """process_fastq_single_end_read_file functions as expected for non-barcoded lane
        """
        actual = process_fastq_single_end_read_file_no_barcode(
            self.fastq1,
            's1',
            min_per_read_length_fraction=0.45)
        actual = list(actual)
        expected = self.fastq1_expected_single_barcode
        self.assertEqual(len(actual), len(expected))
        for i in range(len(expected)):
            np.testing.assert_equal(actual[i], expected[i])

    def test_process_fastq_single_end_read_file_w_defaults_v180(self):
        """process_fastq_single_end_read_file functions as expected w default filters on casava 180 data
        """
        actual = process_fastq_single_end_read_file(self.fastq2,
                                                    self.barcode_fastq2,
                                                    self.barcode_map1,
                                                    min_per_read_length_fraction=0.45)
        actual = list(actual)
        expected = self.fastq2_expected_default
        self.assertEqual(len(actual), len(expected))
        for i in range(len(expected)):
            np.testing.assert_equal(actual[i], expected[i])

    def test_process_fastq_single_end_read_file_handles_log(self):
        """ process_fastq_single_end_read_file generates log when expected
        """
        log = FakeFile()
        list(process_fastq_single_end_read_file(self.fastq1,
                                                self.barcode_fastq1,
                                                self.barcode_map1,
                                                min_per_read_length_fraction=0.45,
                                                log_f=log))
        self.assertTrue(log.s.startswith("Quality filter results"))

    def test_process_fastq_single_end_read_file_handles_histogram(self):
        """ process_fastq_single_end_read_file generates histogram when expected
        """
        histogram = FakeFile()
        list(process_fastq_single_end_read_file(self.fastq1,
                                                self.barcode_fastq1,
                                                self.barcode_map1,
                                                min_per_read_length_fraction=0.45,
                                                histogram_f=histogram))
        self.assertTrue(histogram.s.startswith("Length"))

    def test_check_header_match_pre180(self):
        """check_header_match_pre180 functions as expected with varied input """

        # match w illumina qual string
        self.assertTrue(check_header_match_pre180("@990:2:4:11272:5533#1/1",
                                                  "@990:2:4:11272:5533#1/2"))
        self.assertTrue(check_header_match_pre180("@990:2:4:11272:5533#1/1",
                                                  "@990:2:4:11272:5533#1/3"))
        # qual string differs (this is acceptable)
        self.assertTrue(check_header_match_pre180("@990:2:4:11272:5533#1/1",
                                                  "@990:2:4:11272:5533#0/3"))
        # match wo illumina qual string
        self.assertTrue(check_header_match_pre180("@990:2:4:11272:5533/1",
                                                  "@990:2:4:11272:5533/2"))
        self.assertTrue(check_header_match_pre180("@990:2:4:11272:5533/1",
                                                  "@990:2:4:11272:5533/3"))

        # mismatch w illumina qual string
        self.assertFalse(check_header_match_pre180("@990:2:4:11272:5533#1/1",
                                                   "@990:2:4:11272:5532#1/2"))
        self.assertFalse(check_header_match_pre180("@990:2:4:11272:5533#1/1",
                                                   "@890:2:4:11272:5533#1/2"))
        # mismatch wo illumina qual string
        self.assertFalse(check_header_match_pre180("@990:2:4:11272:5533/1",
                                                   "@990:2:4:11272:5532/2"))
        self.assertFalse(check_header_match_pre180("@990:2:4:11272:5533/1",
                                                   "@890:2:4:11272:5533/2"))

    def test_check_header_match_180_or_later(self):
        """check_header_match_180_or_later functions as expected with varied input """
        # identical
        self.assertTrue(check_header_match_180_or_later(
            "M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0",
            "M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0"))
        # identical except read number
        self.assertTrue(check_header_match_180_or_later(
            "M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0",
            "M00176:17:000000000-A0CNA:1:1:15487:1773 2:N:0:0"))
        # identical except read number
        self.assertTrue(check_header_match_180_or_later(
            "M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0",
            "M00176:17:000000000-A0CNA:1:1:15487:1773 3:N:0:0"))
        # different reads
        self.assertFalse(check_header_match_180_or_later(
            "M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0",
            "M00176:17:000000000-A0CNA:1:1:16427:1774 1:N:0:0"))

    def test_process_fastq_single_end_read_file_toggle_store_unassigned(self):
        """process_fastq_single_end_read_file handles store_unassigned
        """
        fastq_f = [
            "@990:2:4:11272:5533#1/1",
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            "+",
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
            "@990:2:4:11272:5533#1/2",
            "GAAAAAAAAAAT",
            "+",
            "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA': 's1'}
        # empty results when store_unassigned=False
        actual = process_fastq_single_end_read_file(fastq_f,
                                                    barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = []
        self.assertEqual(actual, expected)

        # non-empty results when store_unassigned=True
        actual = process_fastq_single_end_read_file(fastq_f,
                                                    barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=True,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = [(
            'Unassigned_0 990:2:4:11272:5533#1/1 orig_bc=GAAAAAAAAAAT new_bc=GAAAAAAAAAAT bc_diffs=0',
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
           np.array([34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
                     34, 34, 34, 34, 34, 34, 34, 34, 25, 32, 32, 28, 32, 34, 34, 34, 34,
                     34, 34, 34, 34, 34, 34, 34, 34, 34, 32, 34, 34, 34, 34, 33, 34, 32,
                     33, 32, 31, 27, 34, 33, 31, 33, 33, 29, 34, 30, 31, 34,  9, 23, 20,
                     20, 17, 30, 25, 18, 30, 21, 32], dtype=np.int8) ,
            0)]
        np.testing.assert_equal(actual, expected)

    def test_process_fastq_single_end_read_file_toggle_thirteen_base_barcodes(
            self):
        """process_fastq_single_end_read_file handles thriteen base reads of tweleve base barcodes
        """
        fastq_f = [
            "@990:2:4:11272:5533#1/1",
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            "+",
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
            "@990:2:4:11272:5533#1/2",
            "AAAAAAAAAAAAT",
            "+",
            "bbbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA': 's1', 'TAAAAAAAAAAA': 's2'}

        # rev_comp = False
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = [(
            's1_0 990:2:4:11272:5533#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0',
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            np.array([34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 25, 32, 32, 28, 32, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 34, 32, 34, 34, 34, 34, 33, 34, 32,
                      33, 32, 31, 27, 34, 33, 31, 33, 33, 29, 34, 30, 31, 34,  9, 23, 20,
                      20, 17, 30, 25, 18, 30, 21, 32], dtype=np.int8),
            0)]
        np.testing.assert_equal(actual, expected)

    def test_process_fastq_single_end_read_file_toggle_rev_comp(self):
        """process_fastq_single_end_read_file handles rev_comp
        """
        fastq_f = [
            "@990:2:4:11272:5533#1/1",
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            "+",
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
            "@990:2:4:11272:5533#1/2",
            "AAAAAAAAAAAA",
            "+",
            "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA': 's1'}

        # rev_comp = False
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = [(
            's1_0 990:2:4:11272:5533#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0',
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            np.array([34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 25, 32, 32, 28, 32, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 34, 32, 34, 34, 34, 34, 33, 34, 32,
                      33, 32, 31, 27, 34, 33, 31, 33, 33, 29, 34, 30, 31, 34,  9, 23, 20,
                      20, 17, 30, 25, 18, 30, 21, 32], dtype=np.int8),
            0)]
        np.testing.assert_equal(actual, expected)

        # rev_comp = True
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=True,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = [(
            's1_0 990:2:4:11272:5533#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0',
            "GGCTGGCTCCCCTTTCGGGGTTACCTCACCGACTTCGGGTGTTGCCGACTCTCGTGGTGTGACGGGCGGTGTGTGC",
            _ascii_to_phred64("`U^RY^QTTWIb_^b]aa_ab[_`a`babbbb`bbbbbbbbbbbbb`\``Ybbbbbbbbbbbbbbbbbbbbbbbbb"),
            0)]
        np.testing.assert_equal(actual, expected)

    def test_process_fastq_single_end_read_file_error_on_header_mismatch(self):
        """ValueError on barcode/read header mismatch
        """
        fastq_f = [
            "@990:2:4:11272:5533#1/1",
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            "+",
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
            "@990:2:4:11272:5532#1/2",
            "TTTTTTTTTTTT",
            "+",
            "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA': 's1'}
        actual = process_fastq_single_end_read_file(
            fastq_f,
            barcode_fastq_f,
            barcode_to_sample_id,
            store_unassigned=False,
            max_bad_run_length=0,
            phred_quality_threshold=2,
            min_per_read_length_fraction=0.75,
            rev_comp=False,
            rev_comp_barcode=False,
            seq_max_N=0,
            start_seq_id=0)
        self.assertRaises(FastqParseError, list, actual)

    def test_process_fastq_single_end_read_file_toggle_rev_comp_barcode(self):
        """process_fastq_single_end_read_file handles rev_comp_barcode
        """
        fastq_f = [
            "@990:2:4:11272:5533#1/1",
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            "+",
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
            "@990:2:4:11272:5533#1/2",
            "TTTTTTTTTTTT",
            "+",
            "bbbbbbbbbbbb"]
        barcode_to_sample_id = {'AAAAAAAAAAAA': 's1'}
        # empty results when rev_comp_barcode=False
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = []
        self.assertEqual(actual, expected)

        # non-empty results when rev_comp_barcode=True
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=True,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = [(
            's1_0 990:2:4:11272:5533#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0',
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            np.array([34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 25, 32, 32, 28, 32, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 34, 32, 34, 34, 34, 34, 33, 34, 32,
                      33, 32, 31, 27, 34, 33, 31, 33, 33, 29, 34, 30, 31, 34,  9, 23, 20,
                      20, 17, 30, 25, 18, 30, 21, 32], dtype=np.int8),
            0)]
        np.testing.assert_equal(actual, expected)

        # forward orientation no longer matches when rev_comp_barcode=True
        barcode_to_sample_id = {'TTTTTTTTTTTT': 's1'}
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=True,
                                                    seq_max_N=0,
                                                    start_seq_id=0)
        actual = list(actual)
        expected = []
        self.assertEqual(actual, expected)

    def test_process_fastq_single_end_read_file_w_golay_correction(self):
        """process_fastq_single_end_read_file handles golay correction
        """
        fastq_f = [
            "@990:2:4:11272:5533#1/1",
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            "+",
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"]
        barcode_fastq_f = [
            "@990:2:4:11272:5533#1/2",
            "ACAGACCACTCT",
            "+",
            "bbbbbbbbbbbb"]

        barcode_to_sample_id = {'ACAGACCACTCA': 's1'}
        # empty result with single barcode error and golay correction
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0,
                                                    barcode_correction_fn=decode_golay_12,
                                                    max_barcode_errors=1.5)
        actual = list(actual)
        expected = [(
            's1_0 990:2:4:11272:5533#1/1 orig_bc=ACAGACCACTCT new_bc=ACAGACCACTCA bc_diffs=1',
            "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
            np.array([34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 25, 32, 32, 28, 32, 34, 34, 34, 34,
                      34, 34, 34, 34, 34, 34, 34, 34, 34, 32, 34, 34, 34, 34, 33, 34, 32,
                      33, 32, 31, 27, 34, 33, 31, 33, 33, 29, 34, 30, 31, 34,  9, 23, 20,
                      20, 17, 30, 25, 18, 30, 21, 32], dtype=np.int8),
                   0)]
        np.testing.assert_equal(actual, expected)

        # empty result with adjusted max_barcode_errors
        actual = process_fastq_single_end_read_file(fastq_f, barcode_fastq_f,
                                                    barcode_to_sample_id,
                                                    store_unassigned=False,
                                                    max_bad_run_length=0,
                                                    phred_quality_threshold=2,
                                                    min_per_read_length_fraction=0.75,
                                                    rev_comp=False,
                                                    rev_comp_barcode=False,
                                                    seq_max_N=0,
                                                    start_seq_id=0,
                                                    barcode_correction_fn=decode_golay_12,
                                                    max_barcode_errors=0.9)
        actual = list(actual)
        expected = []
        self.assertEqual(actual, expected)

    def test_bad_chars_from_threshold(self):
        """bad_chars_from_threshold selects correct chars as bad
        """
        exp1 = [
            '\t',
            '\n',
            '\r',
            ' ',
            '!',
            '"',
            '#',
            '$',
            '%',
            '&',
            "'",
            '(',
            ')',
            '*',
            '+',
            ',',
            '-',
            '.',
            '/',
            '0',
            '1',
            '2',
            '3',
            '4',
            '5',
            '6',
            '7',
            '8',
            '9',
            ':',
            ';',
            '<',
            '=',
            '>',
            '?',
            '@',
            'A',
            'B']
        exp2 = ['\t',
                '\n',
                '\r',
                ' ',
                '!',
                '"',
                '#',
                '$',
                '%',
                '&',
                "'",
                '(',
                ')',
                '*',
                '+',
                ',',
                '-',
                '.',
                '/',
                '0',
                '1',
                '2',
                '3',
                '4',
                '5',
                '6',
                '7',
                '8',
                '9',
                ':',
                ';',
                '<',
                '=',
                '>',
                '?',
                '@',
                'A',
                'B',
                'C',
                'D',
                'E',
                'F',
                'G',
                'H',
                'I',
                'J',
                'K',
                'L',
                'M',
                'N',
                'O',
                'P',
                'Q',
                'R',
                'S',
                'T',
                'U',
                'V',
                'W',
                'X',
                'Y',
                'Z',
                '[',
                '\\',
                ']',
                '^',
                '_',
                '`',
                'a',
                'b',
                'c',
                'd',
                'e',
                'f',
                'g',
                'h',
                'i',
                'j',
                'k',
                'l',
                'm',
                'n',
                'o',
                'p',
                'q',
                'r',
                's',
                't',
                'u',
                'v',
                'w',
                'x',
                'y',
                'z',
                '{',
                '|',
                '}',
                '~']
        exp3 = [
            '\t',
            '\n',
            '\r',
            ' ',
            '!',
            '"',
            '#',
            '$',
            '%',
            '&',
            "'",
            '(',
            ')',
            '*',
            '+',
            ',',
            '-',
            '.',
            '/',
            '0',
            '1',
            '2',
            '3',
            '4',
            '5',
            '6',
            '7',
            '8',
            '9',
            ':',
            ';',
            '<',
            '=',
            '>',
            '?',
            '@']
        self.assertEqual(bad_chars_from_threshold('B'),
                         {}.fromkeys(exp1))
        self.assertEqual(bad_chars_from_threshold(''), {})
        self.assertEqual(bad_chars_from_threshold('~'),
                         {}.fromkeys(exp2))
        self.assertEqual(bad_chars_from_threshold('@'),
                         {}.fromkeys(exp3))

    def test_quality_filter_sequence_pass(self):
        """quality_filter_sequence functions as expected for good read
        """
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

    def test_quality_filter_illumina_qual(self):
        """quality_filter_sequence functions as expected with bad illumina qual digit
        """
        # header with no qual data passes
        header = "990:2:4:11271:5323/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=0.75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # header with no qual data passes
        header = "990:2:4:11271:5323/0"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # header with no qual data passes (old barcode in header format)
        header = "HWI-6X_9267:1:1:4:1699#ACCACCC/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # bad qual fails filter
        header = "@HWI-ST753_50:6:1101:1138:1965#0/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (3,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # bad qual passes filter if filter turned off
        header = "@HWI-ST753_50:6:1101:1138:1965#0/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=False)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # good qual passes filter
        header = "@HWI-ST753_50:6:1101:1138:1965#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

    def test_quality_filter_sequence_fail_w_B(self):
        """quality_filter_sequence handles bad qual score as expected
        """

        # early 'B' in sequence causes truncation and too short of a read
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            _ascii_to_phred64("bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        np.testing.assert_equal(
            actual,
            (1,
             "GCACTCACCGCCCGTCAC",
             _ascii_to_phred64("bbbbbbbbbbbbbbbbbb")))

        # increasing max_bad_run_length rescues read
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=1,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # changing threshold rescues read
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=1,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        self.assertEqual(actual, (0,
                                  "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                                  "bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"))

        # changing min_per_read_length_fraction rescues read
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            _ascii_to_phred64("bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=5,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        np.testing.assert_equal(
            actual,
            (0,
             "GCACTCACCGCCCGTCAC",
             _ascii_to_phred64("bbbbbbbbbbbbbbbbbb")))

    def test_quality_filter_sequence_fail_w_N(self):
        """quality_filter_sequence handles N as expected
        """

        # 'N' in sequence causes failure
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)
        expected = (2,
                    "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                    "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
        self.assertEqual(actual, expected)

        # increasing max N rescues sequence
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC"
        quality =  \
            "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=75,
                                         seq_max_N=1,
                                         filter_bad_illumina_qual_digit=True)

        expected = (0,
                    "GCACTCACCGCCCGTCACACCACGAAAGTNGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
                    "bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
        self.assertEqual(actual, expected)

        # truncation of N rescues sequence (sequence is truncated when
        # the quality hits B, and the truncated sequence is above the
        # length threshold and no longer contains an N)
        header = "990:2:4:11271:5323#1/1"
        sequence = \
            "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTN"
        quality =  \
            _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^B`")
        actual = quality_filter_sequence(header,
                                         sequence,
                                         quality,
                                         max_bad_run_length=0,
                                         phred_quality_threshold=2,
                                         min_per_read_length=50,
                                         seq_max_N=0,
                                         filter_bad_illumina_qual_digit=True)

        expected = (0,
                    "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTG",
                    _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^"))
        np.testing.assert_equal(actual, expected)

barcode_map1 = {'AAAAAAAAAAAA': 's1',
                'AAAAAAAAAAAC': 's2',
                'AAAAAAAAAAAG': 's3',
                'AAAAAAAAAAAT': 's4', }

fastq1 = """@990:2:4:11271:5323#1/1
GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC
+
bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`
@990:2:4:11271:5323#1/1
GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG
+
bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT
@990:2:4:11272:9538#1/1
GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG
+
b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa
@990:2:4:11272:9538#1/1
GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC
+
bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O
@990:2:4:11272:7447#1/1
GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:7447#1/1
GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC
+
b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:19991#1/1
GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:19991#1/1
GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA
+
bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB
@990:2:4:11272:4315#1/1
GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG
+
bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:4315#1/1
GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG
+
``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:5533#1/1
GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC
+
``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@990:2:4:11272:5533#0/1
GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC
+
bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB
"""

barcode_fastq1 = """@990:2:4:11271:5323#1/2
AAAAAAAAAAAA
+
bbbbbbbbbbbb
@990:2:4:11271:5323#1/2
AAAAAAAAAAAC
+
bbcbbbbbbbbb
@990:2:4:11272:9538#1/2
AAAAAAAAAAAA
+
b_bbbbbbbbbb
@990:2:4:11272:9538#1/2
AAAAAAAAAAAT
+
bb^bbbbbbbbb
@990:2:4:11272:7447#1/2
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@990:2:4:11272:7447#1/2
AAAAAAAAAAAA
+
b`bbbbbbbbbb
@990:2:4:11272:19991#1/2
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@990:2:4:11272:19991#1/2
AAAAAAAAAAAC
+
bbbbbbbbbbbb
@990:2:4:11272:4315#1/2
AAAAAAAAAAAT
+
bbbb_bbbbbbb
@990:2:4:11272:4315#1/2
AAAAAAAAAAAT
+
``Q``````_``
@990:2:4:11272:5533#1/2
GAAAAAAAAAAT
+
``Q``````_``
@990:2:4:11272:5533#0/2
AAAAAAAAAAAT
+
bbbbbbbbbbbb
"""

fastq2 = """@M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0
GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC
+
CCCCCCCCCC###############:AA=ACCCCCCCCCCCCCACCCCBCABA@<CB@BB>C?@C*8552?:3?6A
@M00176:17:000000000-A0CNA:1:1:17088:1773 1:N:0:0
GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG
+
CCDCCCCCCCCCCCCCCCCCCCCCCCCCC@CCCCCCCCBCB@C?C:@ABB?C1CACCCC):(:;5CC?@BC<?CB5
@M00176:17:000000000-A0CNA:1:1:16738:1773 1:N:0:0
GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG
+
C@CCCCCCCCCCCCCCCCCCCCCCCCCCBCBB?BA<CCCCACCCC5CCBCC>C><@BAB>BDBBBDCBDB@B?ABB
@M00176:17:000000000-A0CNA:1:1:12561:1773 1:N:0:0
GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC
+
CC?CCC####CCCCCCCCCCCCCCCCBCCCCAACCCA@@CCCCCC*8399A3AA=A=:=?@@CB?B<4BBB@>0>0
@M00176:17:000000000-A0CNA:1:1:14596:1773 1:N:0:0
GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG
+
CACCCCCCCCCCCCCCCA?CCCCC:CCCCC=@@@A@CCBC?BBB6?=A############################
@M00176:17:000000000-A0CNA:1:1:12515:1774 1:N:0:0
GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC
+
CACCCCCCCCCCCCCCCA?CCCCC:CCCCC=@@@A@CCBC?BBB6?=A############################
@M00176:17:000000000-A0CNA:1:1:17491:1774 1:N:0:0
GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG
+
CCCCCCCCCCCCCCCCCCCCC9CCC@CCCBCCCAB;<6>=05:97A5C############################
@M00176:17:000000000-A0CNA:1:1:16427:1774 1:N:0:0
GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA
+
CCCCCCCCCCCCCCCCCCCCBACCCCCCCCCCABCC@BBDCCCCC>@@@>=<=??<B0D#################
@M00176:17:000000000-A0CNA:1:1:13372:1775 1:N:0:0
GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG
+
CCCC@CCCCCCCCCCAAA2AAA######################################################
@M00176:17:000000000-A0CNA:1:1:14806:1775 1:N:0:0
GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG
+
AA2AAAAAA@AA####AAAA,>>B####################################################
@M00176:17:000000000-A0CNA:1:1:13533:1775 1:N:0:0
GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC
+
AA2AAAAAA@AAAAAAAAAA,>>B####################################################
@M00176:17:000000000-A0CNA:1:1:18209:1775 1:N:0:0
GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC
+
CCCCCCCCCCCCCCCCCCCCC9CCC@CCCBCCCAB;<6>=05:97A5C############################
"""

barcode_fastq2 = """@M00176:17:000000000-A0CNA:1:1:15487:1773 2:N:0:0
AAAAAAAAAAAA
+
AAAAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:17088:1773 2:N:0:0
AAAAAAAAAAAC
+
AABAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:16738:1773 2:N:0:0
AAAAAAAAAAAA
+
A>AAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:12561:1773 2:N:0:0
AAAAAAAAAAAT
+
AA:AAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:14596:1773 2:N:0:0
AAAAAAAAAAAA
+
A?AAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:12515:1774 2:N:0:0
AAAAAAAAAAAA
+
A?AAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:17491:1774 2:N:0:0
AAAAAAAAAAAC
+
AAAAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:16427:1774 2:N:0:0
AAAAAAAAAAAC
+
AAAAAAAAAAAA
@M00176:17:000000000-A0CNA:1:1:13372:1775 2:N:0:0
AAAAAAAAAAAT
+
AAAA>AAAAAAA
@M00176:17:000000000-A0CNA:1:1:14806:1775 2:N:0:0
AAAAAAAAAAAT
+
>>1>>>>>>;>>
@M00176:17:000000000-A0CNA:1:1:13533:1775 2:N:0:0
GAAAAAAAAAAT
+
>>1>>>>>>;>>
@M00176:17:000000000-A0CNA:1:1:18209:1775 2:N:0:0
AAAAAAAAAAAT
+
AAAAAAAAAAAA
"""


fastq1_expected_no_qual_unassigned = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s2_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     _ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     _ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s4_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     _ascii_to_phred64("bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     5),
    ("s2_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     6),
    ("s2_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB"),
     7),
    ("s4_8 990:2:4:11272:4315#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG",
     _ascii_to_phred64("bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     8),
    ("s4_9 990:2:4:11272:4315#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG",
     _ascii_to_phred64("``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     9),
    ("Unassigned_10 990:2:4:11272:5533#1/1 orig_bc=GAAAAAAAAAAT new_bc=GAAAAAAAAAAT bc_diffs=0",
     "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
     _ascii_to_phred64("``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     10),
    ("s4_11 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     11)]

fastq1_expected_default = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s2_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     _ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     _ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s4_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     _ascii_to_phred64("bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     5),
    ("s2_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"),
     6),
    ("s2_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc"),
     7),
    ("s4_8 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"), 8)]

fastq1_expected_single_barcode = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s1_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     _ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     _ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s1_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     _ascii_to_phred64("bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     5),
    ("s1_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"),
     6),
    ("s1_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc"),
     7),
    ("s1_8 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"), 8)]

fastq2_expected_default = [
    ("s1_0 M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     _ascii_to_phred64("bbbbbbbbbbBBBBBBBBBBBBBBBY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s2_1 M00176:17:000000000-A0CNA:1:1:17088:1773 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     _ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 M00176:17:000000000-A0CNA:1:1:16738:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     _ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s4_3 M00176:17:000000000-A0CNA:1:1:12561:1773 1:N:0:0 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     _ascii_to_phred64("bb^bbbBBBBbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 M00176:17:000000000-A0CNA:1:1:14596:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     4),
    ("s1_5 M00176:17:000000000-A0CNA:1:1:12515:1774 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     _ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     5),
    ("s2_6 M00176:17:000000000-A0CNA:1:1:17491:1774 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"),
     6),
    ("s2_7 M00176:17:000000000-A0CNA:1:1:16427:1774 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc"),
     7),
    ("s4_8 M00176:17:000000000-A0CNA:1:1:18209:1775 1:N:0:0 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     _ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"), 8)]


if __name__ == "__main__":
    main()
