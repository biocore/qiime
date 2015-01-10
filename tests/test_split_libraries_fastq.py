#!/usr/bin/env python
# File created on 05 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.0-rc2"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import numpy as np

from unittest import TestCase, main
from tempfile import mkdtemp, NamedTemporaryFile
from shutil import rmtree

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
    process_fastq_single_end_read_file_no_barcode,
    extract_reads_from_interleaved
)
from qiime.golay import decode_golay_12

import skbio.parse.sequences
from skbio.parse.sequences.fastq import ascii_to_phred64, ascii_to_phred33


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

        # vars for test_create_forward_and_reverse_fp
        self.temp_dir_path = mkdtemp()
        self.create_forward_and_reverse = NamedTemporaryFile(
            prefix='create_forward_and_reverse_fp_',
            suffix='.fastq',
            dir=self.temp_dir_path,
            delete=False)
        self.create_forward_and_reverse_fp = self.create_forward_and_reverse.name
        self.create_forward_and_reverse.write(forward_reads)
        self.create_forward_and_reverse.write(reverse_reads)
        self.create_forward_and_reverse.close()

    def tearDown(self):
        """Remove all temp files"""
        rmtree(self.temp_dir_path)

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

    def test_process_fastq_single_end_read_file_invalid_phred_offset(self):
        # passing phred_offset that isn't 33 or 64 raises error
        with self.assertRaises(ValueError):
            list(process_fastq_single_end_read_file(
                self.fastq1,self.barcode_fastq1, self.barcode_map1,
                store_unassigned=True, max_bad_run_length=1000,
                phred_quality_threshold=None, min_per_read_length_fraction=0.,
                rev_comp=False, rev_comp_barcode=False, seq_max_N=1000,
                start_seq_id=0, filter_bad_illumina_qual_digit=False,
                phred_offset=42))

        # passing wrong phred_offset for illumina1.8+ data raises error
        with self.assertRaises(skbio.parse.sequences.FastqParseError):
            list(process_fastq_single_end_read_file(
                self.fastq2, self.barcode_fastq2, self.barcode_map1,
                min_per_read_length_fraction=0.45, phred_offset=64))

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
        # test autodetection of phred_offset (phred_offset=None) and
        # phred_offset=33 gives same results
        for phred_offset in None, 33:
            actual = process_fastq_single_end_read_file(
                self.fastq2, self.barcode_fastq2, self.barcode_map1,
                min_per_read_length_fraction=0.45, phred_offset=phred_offset)

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
            ascii_to_phred64("`U^RY^QTTWIb_^b]aa_ab[_`a`babbbb`bbbbbbbbbbbbb`\``Ybbbbbbbbbbbbbbbbbbbbbbbbb"),
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
            ascii_to_phred64("bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
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
             ascii_to_phred64("bbbbbbbbbbbbbbbbbb")))

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
            ascii_to_phred64("bbbbbbbbbbbbbbbbbbBbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`")
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
             ascii_to_phred64("bbbbbbbbbbbbbbbbbb")))

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
            ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^B`")
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
                    ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^"))
        np.testing.assert_equal(actual, expected)

    def test_create_forward_and_reverse_fp(self):
        """ perform different tests for extract_reads_from_interleaved """

        # regular processing
        extract_reads_from_interleaved(
            self.create_forward_and_reverse_fp, '1:N:0', '2:N:0',
            self.temp_dir_path)
        forward = open(self.temp_dir_path + '/forward_reads.fastq', 'U').read()
        reverse = open(self.temp_dir_path + '/reverse_reads.fastq', 'U').read()

        self.assertEqual(forward, forward_reads)
        self.assertEqual(reverse, reverse_reads)

        # should raise an error due to no matching id
        with self.assertRaises(ValueError):
            extract_reads_from_interleaved(
                self.create_forward_and_reverse_fp, '1N', '2N',
                self.temp_dir_path)


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
CCCCCCCCCCAAAAAAAAAAAAAAA:AA=ACCCCCCCCCCCCCACCCCBCABA@<CB@BB>C?@C*8552?:3?6A
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
CC?CCCAAAACCCCCCCCCCCCCCCCBCCCCAACCCA@@CCCCCC*8399A3AA=A=:=?@@CB?B<4BBB@>0>0
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
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s2_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s4_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     ascii_to_phred64("bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGGCAGTCTAACCGCAAGGAGGACGCTGTCG",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACAGCGTCCTCCTTGCTGTTAGACTTCCGGC",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`BBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     5),
    ("s2_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGGTGGGGCAACCGGCTGTCCCTTTTAGCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     6),
    ("s2_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGGCTGGCCCTTTCCACCCA",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOcBBBBBBBBBBBBBBBBB"),
     7),
    ("s4_8 990:2:4:11272:4315#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GTACTCACCGCCCGTCACGCCATGGGAGTTGGGCTTACCTGAAGCCCGCGAGCTAACCGGAAAGGGGGGGATGTGG",
     ascii_to_phred64("bbbb_bbbbbbbbbb```Q```BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     8),
    ("s4_9 990:2:4:11272:4315#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCCCGTCGCTCGGCGTACCTTCGACCGCTGCCTCCTTTTGGTTATATCTCCGGG",
     ascii_to_phred64("``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     9),
    ("Unassigned_10 990:2:4:11272:5533#1/1 orig_bc=GAAAAAAAAAAT new_bc=GAAAAAAAAAAT bc_diffs=0",
     "GCACACACCGCCCGTCACACCACGAGAGTCGGCAACACCCGAAGTCGGTGAGGTAACCCCGAAAGGGGAGCCAGCC",
     ascii_to_phred64("``Q``````_``````````K]]aBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     10),
    ("s4_11 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCGGCTGGCTCCCCTTTCGGGGGTACCTCAC",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`TbBBBBBBBBBBBBBBBBBBBBBBBBBBBB"),
     11)]

fastq1_expected_default = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s2_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s4_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     ascii_to_phred64("bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     5),
    ("s2_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"),
     6),
    ("s2_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc"),
     7),
    ("s4_8 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"), 8)]

fastq1_expected_single_barcode = [
    ("s1_0 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbbbbbY``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s1_1 990:2:4:11271:5323#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s1_3 990:2:4:11272:9538#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     ascii_to_phred64("bb^bbbbbbbbbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     4),
    ("s1_5 990:2:4:11272:7447#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     5),
    ("s1_6 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"),
     6),
    ("s1_7 990:2:4:11272:19991#1/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc"),
     7),
    ("s1_8 990:2:4:11272:5533#0/1 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"), 8)]

fastq2_expected_default = [
    ("s1_0 M00176:17:000000000-A0CNA:1:1:15487:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACTCACCGCCCGTCACACCACGAAAGTTGGTAACACCCGAAGCCGGTGAGATAACCTTTTAGGAGTCAGCTGTC",
     ascii_to_phred64("bbbbbbbbbb```````````````Y``\`bbbbbbbbbbbbb`bbbbab`a`_[ba_aa]b^_bIWTTQ^YR^U`"),
     0),
    ("s2_1 M00176:17:000000000-A0CNA:1:1:17088:1773 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGTTACCTTGTTACGACTTCACCCCAATCATCGGCCCCACCTTAGACAGCTGACTCCTAAAAGGTTATCTCACCGG",
     ascii_to_phred64("bbcbbbbbbbbbbbbbbbbbbbbbbbbbb_bbbbbbbbaba_b^bY_`aa^bPb`bbbbHYGYZTbb^_ab[^baT"),
     1),
    ("s1_2 M00176:17:000000000-A0CNA:1:1:16738:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGAGGTACCCGAAGCCGGTAGTCTAACCGCAAGGAGGACGCTGTCG",
     ascii_to_phred64("b_bbbbbbbbbbbbbbbbbbbbbbbbbbabaa^a`[bbbb`bbbbTbbabb]b][_a`a]acaaacbaca_a^`aa"),
     2),
    ("s4_3 M00176:17:000000000-A0CNA:1:1:12561:1773 1:N:0:0 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGCTACCTTGTTACGACTTCACCCTCCTCACTAAACGTACCTTCGACAGCGTCCTCCTTGCGGTTAGACTACCGGC",
     ascii_to_phred64("bb^bbb````bbbbbbbbbbbbbbbbabbbb``bbb`__bbbbbbIWRXX`R``\`\Y\^__ba^a[Saaa_]O]O"),
     3),
    ("s1_4 M00176:17:000000000-A0CNA:1:1:14596:1773 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GCACACACCGCCCGTCACACCATCCGAGTTGGGGGTACCCGAAGCCGG",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     4),
    ("s1_5 M00176:17:000000000-A0CNA:1:1:12515:1774 1:N:0:0 orig_bc=AAAAAAAAAAAA new_bc=AAAAAAAAAAAA bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCTCCTCACTCATCGTACCCTCGACA",
     ascii_to_phred64("b`bbbbbbbbbbbbbbb`^bbbbbYbbbbb\___`_bbab^aaaU^\`"),
     5),
    ("s2_6 M00176:17:000000000-A0CNA:1:1:17491:1774 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GCACTCACCGCCCGTCACGCCACGGAAGCCGGCTGCACCTGAAGCCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"),
     6),
    ("s2_7 M00176:17:000000000-A0CNA:1:1:16427:1774 1:N:0:0 orig_bc=AAAAAAAAAAAC new_bc=AAAAAAAAAAAC bc_diffs=0",
     "GGCTACCTTGTTACGACTTCGCCCCAGTCACCGACCACACCCTCGACGGCTGCCTCCGG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbba`bbbbbbbbbb`abb_aacbbbbb]___]\[\^^[aOc"),
     7),
    ("s4_8 M00176:17:000000000-A0CNA:1:1:18209:1775 1:N:0:0 orig_bc=AAAAAAAAAAAT new_bc=AAAAAAAAAAAT bc_diffs=0",
     "GGATACCTTGTTACGACTTCACCCCAATCATCGACCCCACCTTCGGCG",
     ascii_to_phred64("bbbbbbbbbbbbbbbbbbbbbXbbb_bbbabbb`aZ[U]\OTYXV`Tb"), 8)]


forward_reads = """@MISEQ03:64:000000000-A2H3D:1:1101:14358:1530 1:N:0:TCCACAGGAGT
TNCAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGAAACTGACACTGAGGGGCGAAAGCGGGGGGGGCAAACG
+
?#5<????DDDDDDDDEEEEFFHHHHHHHHHHHHHHDCCHHFGDEHEH>CCE5AEEHHHHHHHHHHHHHHHHHFFFFHHHHHHEEADEEEEEEEEEEEEEEEEEEEEEEE?BEEEEEEEEEEEAEEEE0?A:?EE)8;)0ACEEECECCECAACEE?>)8CCC?CCA8?88ACC*A*::A??:0?C?.?0:?8884>'.''..'0?8C?C**0:0::?ECEE?############################
@MISEQ03:64:000000000-A2H3D:1:1101:14206:1564 1:N:0:TCCACAGGAGT
TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGTCAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGTTTGAAACTACAAGGCTAGAGTGTAGCAGAGGGGGGTAGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGGGGAGGAATACCAATGGCGAAGGCAGCCCCCGGGGTTAACACTGACGCCAAGGCACGAAAGCGGGGGGGGCAAACG
+
?????BB?DDDDDD@DDCEEFFH>EEFHHHHHHGHHHCEEFFDC5EECCCCCCDECEHF;?DFDDFHDDBBDF?CFDCCFEA@@::;EEEEEEEECBA,BBE?:>AA?CA*:**0:??A:8*:*0*0**0*:?CE?DD'..0????:*:?*0?EC*'.)4.?A***00)'.00*0*08)8??8*0:CEE*0:082.4;**?AEAA?#############################################
@MISEQ03:64:000000000-A2H3D:1:1101:14943:1619 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCGACCACCTGGACTGAAACTGACACTGAGGGGCGAAAGCGGGGGGGGCAAAAG
+
?AAAABBADDDEDEDDGGGGGGIIHHHIIIIIIIHHHCCHHFFEFHHHDCDH5CFHIIHIIIIHHHHHHHHHHHHHHHHHHHGGGEGGGGDDEGEGGGGGGGGGGGGGGEEEGCCGGGGGGEGCEEEECE?ECGE.84.8CEEGGGE:CCCC0:?C<8.48CC:C??.8.8?C:*:0:*9)??CCEE**)0'''42<2C8'8..8801**0*.1*1?:?EEEC?###########################
@MISEQ03:64:000000000-A2H3D:1:1101:15764:1660 1:N:0:TCCACAGGAGT
TACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGATTGTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGGCGATCTTGAGTCCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGGAACTGACGCTGAGGGGCGAAAGCTGGGGGAGCAAACG
+
???????@DDDDDDBDFEEFEFHEHHHHHHHHHHHHHEHHHHFEHHHHAEFHGEHAHHHHHHHHHHHHHHH=@FEFEEFEFEDAEEEFFE=CEBCFFFCECEFEFFFCEEEFFCD>>FFFEFF*?EED;?8AEE08*A*1)?E::???;>2?*01::A?EEEFEEE?:C.8:CE?:?8EE8AECEFE?C0::8'488DE>882)*1?A*8A########################################
@MISEQ03:64:000000000-A2H3D:1:1101:15211:1752 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGGGGTGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATGGGAAGGAACACCAGGGGCGAAGGCGACCACCTGGACTGATACTGACACTGGGGTGGGAAAGGGGGGGGAGGAAAAG
+
?????<B?DBBDDDBACCEEFFHFHHHHHHHHHHHHH5>EFFFEAACEC7>E5AFEFHHHHHHF=GHFGHFHHHHFHFHH;CED8@DDDE=4@EEEEECE=CECECEECCBB34,=CAB,40:?EEEE:?AAAE8'.4'..8*:AEEECCCA::A################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:15201:1774 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGCGAAAGCGGGGGGGGCAAACG
+
?????BB?DDDDDDBDEEEEFFHFHHHHHHHHHHHFH>CEHDGDDCDE5CCEACFHHHHHHHHFFHHHHHHHHFHHFHHHHHHDEBFEEE@DEEEEEEEEEEEEEEBBCBECEEEEEEEEEEEEEEE?ACCEEEA)84)0.?EEE:AEACA?0?CEDD'.4?A:ACA)0'80:A:?*0*0)8CEEEEE?*0*)88888A'.5;2A)*0000*8:*0:?CEEEE?A*?A#######################
@MISEQ03:64:000000000-A2H3D:1:1101:15976:1791 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATGGGAAGGAACACCGGGGGGGAGGGGGGCTCTCGGGTCCTTTTCGGCGGCTGGGGGCGGAAGGCAGGGGGGGCAACCG
+
?????BB?DDDDDDDDEEEEFFIFHHHHHHIIIHIFHCCHF@F>CECHCDDECCFEADEHHHHHHHHFGHHHHHHFHHHHHHF8:<DEEEADEEFFFFFFABEFFEFFECBCEEFEFFFFEACEEEEE*10*1??.08.888AEF?EEEC1:1:??>>'88AC?::?AA##################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:16031:1840 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCGACCACCTGGACGGATACTGACACTGAGGGGCGAAAGGGTGGGGAGAAAAAG
+
?????BB?DDDDDDDDGFEGGGIHHHHIIIIIIIHFE>CFFFFDCHCH>>CE-5EEIIHHHIHHHHHHHHHHGHFDFHFHEHGBEEEEGGEDGGGGEGGEGGGGGCEGCCEEGGG><CEECCGCEEEG?C:1?EG.8<.88CCCEEGE?C?C*:1:<>'.8?8:C:?00.0?:?*1::*9CC?EEEG*?##############################################################
@MISEQ03:64:000000000-A2H3D:1:1101:12964:1852 1:N:0:TCCACAGGAGT
TNCAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGCAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGAAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCAACCACCGGGACTGAAACTGAACCGGAGGGGGGAAAGCGGGGGGGGAAAACG
+
?#55????DDDDDDDDEEEEFFHHBHHHFFGHHFHDC+5>C/?E7DCHCCCD5AECFHHHFHHHHHHHHHFFFFFHFFDFEFF5)@=DEFDEFEEFF;AEAABC,4BECCC=B,5?C8?CC?CC*:?E010:?EA)0.)08C?A:?A########################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17245:1906 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAGGGAACACAAGGGGCGAAGGCGACCACCGGGACGGAAACTGCAACTGGGGGGGGAAAGCGGGGGGGGAAACAG
+
AAA??BB?DDDDDDDDGGEGGGIHGHIIIIIIHF?CFCCDFFFDCHEHC>DH5AFEHIHHHHHHHHHHHHHHFFFFFHHHHHGDBEEGGGGGGG@EGEGGGGGGGCGEGACC>EGEGGGGC:C0CEEG:0::CEE)88)08?:CECCE:C*10*104A.4CE:*:?C8)'8CC##############################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:18171:1983 1:N:0:TCCACAGGAGT
GACGTAGGAGGCGAGCGTTGTCCGGATTCATTGGGCGTAAAGAGCGCGCAGGCGGCTTGGTAAGTCGGATGTGAAATCCCGAGGCTCAACCTCGGGTCTGCATCCGATACTGCCCGGCTAGAGGTAGGTAGGGGAGATCGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGGGGGCGAAGGCGGATCTCTGGGCCTTCCCTGACGCTCAGGCGCGAAAGCGGGGGGGGCGAACG
+
??????B?DDDDDDDDFFEFFFIHFEEEHHIHHHFHHEHHFGFFFHCEHEHCDECCEFFE4DDFDBEEEEEFFFFEEFFCE8B>BEFEEFFCEFE>8>EFFE*A?A?ADDAAEE8E>;>EA:??1*:?111?C<88AA08?ACECF:*:?*08:0:8<.4?EE*A:))'..0*01*?:08?A*?CA#################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:14225:2009 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGAAACGGACACTGAGGGGCGAAAGCGGGGGGGGCAAACG
+
?????BB?DDDDDDBDEEEEFFHHHHIIIIHHIIIIHHEHIFGEHHHHCCEHAEFHIIHIIIIHHHHHHHHHHFHHHHHHHHFFFEFFFFFEFFFFFFEFFFFFFEFFFEFCACEFFFFFF:C?CEEE*?:AAEE88;088?AEFCEAEECEEEFE>?).?ECCEEE8?4AFFE0?*0088ACFFFAAC*0?C888>>CD?D;8CE*:*:A?CF*::)0?DD?:::?########################
@MISEQ03:64:000000000-A2H3D:1:1101:16656:2052 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCGACCACCGGGACTGAAACTGACACTGAGGGGGGAAAGCGGGGGGGGAAAACA
+
?????BB?BDDDDDDDGFFEGGIIHHIIIIIHHHHIHCCFFDEEEHEHFFEH5AFHHIHIHIHGGHHHHHHHFHHFHHHHHHGEG@EGEGGEGGGGCEGGEGGGGEGGACECGGGGGGGGEGGCCEGG?CCCEGC088)0.?EGG?EC*::C*:??<8.48?C:?C808.8CEE#############################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:18209:2060 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTCAGCAAGTCAGGGGTGAAATCCCGGGGCTCAACCCCGGAACTGCCCTTGAAACTGCTAAGCTAGAATCTTGGAGAGGCGAGTGGAATTCCGAGTGTAGAGGGGAAATTCGTAGATATTCGGAAGAACACCAGGGGCGAAGGCGACCCCCTGGACAAGCATTGACGCTGAGGGGGGAAAGCGGGGGGGGCAAAAG
+
?????BB?BDDDDDDDECEEFFHHHHAHFHHHHHHHHCCHHH=@DEEHFHFCGHHB)?ECGHHH?DHHHHHCCCFFHHHFEEEEEEEEEEEEEB)>EDACEECEECEEECEE:*0A:AEAECA:0::ACE??E?.8'4.88?EC*00:08).0:*00?)..8AAAAA*0)0::?::?0A8)?C:?A#################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:13824:2095 1:N:0:TCCACAGGAGT
TACGTAGGGGGCTAGCGTTGTCCGGAATCATTGGGCGTAAAGCGCGTGTAGGCGGCCCGGTAAGTCCGCTGTGAAAGTCGGGGGCTCAACCCTCGAAAGCCGGGGGATACTGTCGGGCTAGAGTACGGAAGAGGCGAGTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCCATTGCGAAGGCAGCTCGCTGGGACGTTACTGAGGCTGAGACCGGAAAGGGGGGGGGGCAAAAG
+
??A??BBADDDDDDBDFEEFFFHHHHHFHHHIHHFHHCCHHFHDCCDEEHHFIHAHHHHH@EFFDFFEBDEDEFFECBBEEEED?28CCFFECE;EF8?ECD;'488?EEFCE:A>>?>EECEE::A8E8.8?8).'.'08AEE*?:*::*001:?<D.'8??*:*))'''01***10*088CEEEEA8C#############################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17839:2106 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACAACCGGGACGGAAACTGACACTGAGGGGCGAAAGCGGGGGGGGCAAAAG
+
AAA?ABB?DDDDDDDEGFEGGFHIHHIIIIIIIIDFH5CFHHGHEH=DC>CE5AEEHFHIHIFHHHHHHHHHFHHFHHHHHHGGGGGEEGGGGGDEGGGGGGGGGGGGGCE>AEGEGGGGEEECEGEE1:??CEC08>.88CEEECG*:C?CC:?0.4.4CE?CECC?)4?CC:*11?:?)CCEGG).9*1:?8<2<<C####################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17831:2124 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGGGGGGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGAAACTGACACTGAGGGGCGAAAGCGGGGGGAGCAAACG
+
AAAAABB?DDDDDDDDGFGFGGIIHHIIIIIHIIDFH>CFHHGDCFDH>CDHAEFEHIEFFGGHHHHHHHFH=CFFHHHHEHG8DEEGEGGGGGDEEEEGEEGGGCGGEEECCACCEGGGCEE::?CE0?CCEGE'.<'..4CEGEGGEEEE*::C>20>?C?*1:C..'8:??*:*?*0)??9CEG8?*1*8'4.44?58<28?C#############################################
@MISEQ03:64:000000000-A2H3D:1:1101:12555:2129 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACGAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACGGAAACTGACACTGAGGTGCGAAAGCGTGGGGACCAACCG
+
????ABBADDDDDEDDGGGGGGIIHHIIIIIHIIHHHFFHHHFHHHHH>CDHAEFHFGHHFHHHHHHHHHFHHFHFFHHHHHGBEEAD;DGGGGEGGGCGCEEEGEGGGCE>>>CEDGDGE:C:CGGG:?C??EE08<)0?ECEGEGCCECEEGGGGG08CECE?CE8)4?CC:*:*:0989*9CEC8C*:?C'842.8'.4.2?E9?*:?'.8).::::?CC:*110*0C8C<8??C#############
@MISEQ03:64:000000000-A2H3D:1:1101:13627:2170 1:N:0:TCCACAGGAGT
GACAGAGGGTGCAAACGTTGTTCGGAATTACTGGGCATAAAGAGCACGTAGGTGGTCTGCTAAGTCACGTGTGAAATCCCCCGGCTCAACCGGGGAATTGCGCGTGATACTGGCCGGCTCGAGGTGGGTAGGGGGGAGCGGAACTCCAGGGGGAGCGGGGAAATGCGTAGATATCTGGAGGAACACCGGGGGCGAAAGCGGCTCACGGGACCCAATCTGACACTGAGGGGCGAAAGCTAGGGTGGCAAACG
+
?????BB?DDDDDDDDEFFFFFHHHHHIHIIHIIFHCEHIIHBFHIHHAAFH5CF@FHHHGHIIGHHHHFHIHIIIHIIIHHHHHHHHHHFHHHFFEFEFEDBE<>BBEEFECECE?D'..2AD)8A>40?AED''''.4<D>>AC**1?).2'888D'5<EACEEEAEDEFEE:*??*08A?AAC)58'4>2<>D8A:A82'.*:*.'?>E)AA####################################
@MISEQ03:64:000000000-A2H3D:1:1101:11781:2223 1:N:0:TCCACAGGAGT
TACGTAGGGCGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTCGCGTCTGCCGTGAAAGTCCGGGGCTCAACTCCGGATCTGCGGTGGGTACGGGCAGACTAGAGTGATGTAGGGGAGATTGGAATTCCTGGTGTAGCGGGGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCATTAACTGACGCTGAGGAGCGAAAGCAGGGGGGGCGAACG
+
???A?BB?DDDDDDDDEEECFFHHHHHIHHHIIIHHHECEHFCFHGGH>CFEFEHHHHHFFDFHCDEFFHHEBFFFF?BBEEEEEEEFFFBEEEEAEDEFEDD.8A8.ACEDDD;AEFFFFEF:*1:?ACCFFD8<AE?EFFFF:EEEEFFFA:CEDD'.8??CEF?ADDDFF:C:?::?AEEFFFA>8'08:2448DE?E?8:*:*1A***0*:AA*?AEEEEE?#########################
@MISEQ03:64:000000000-A2H3D:1:1101:17996:2254 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCGACCACCGGGCCTGAAACTGACACTGAGGGGGGAAAGCGGGGGGGGAAAACG
+
?????BB?DDDDDDDDGGGGGGIHHHHIIIIHHHFFH>CHFHGHHHEHCCCE5AFEHIHHHHHHHHHHHHHHHHHHHHHHHHGGEEGGEGEGGGEGEGGGCGGGGGGGECGEECGAECGGEEEC**CE?C::CCC.8<)08?CCC:CCCEC?CC?:8>'4<?.1C:8082CCGG*:*:0C8?EC*0C89.?############################################################
@MISEQ03:64:000000000-A2H3D:1:1101:13712:2276 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGGGGGAAAGCGGGGGGGGAAAACG
+
?????BBADDDDDDDDGGEGGGIHHHHIIIIIIIIIHCCHHDCECHEHCDEH-AFEHIHIHHIHHHHHHHHHHHFFFHHHHHGEGEDDEEDDDGGGGGEGGGGGEEEGEEGEGGGGGGGCEGEGCEGG:C::CEE)88.8?EGGG:C?:?:C??:*52'.888:CEE).2CCGE*C??:C.?EGGGGC9*118>>.4>C''.8<CC*?*:**00*01?:CEGCC###########################
@MISEQ03:64:000000000-A2H3D:1:1101:15819:2301 1:N:0:TCCACAGGAGT
TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTGAGTTAAGTCTGCTGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGGGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCGACCACCTGGACTGAAACTGACACTGAGGGGCGAAAGCGGGGGGGGCAAACG
+
?????BB?DEDDDDDDFDEEFFHHHEFHHIIIIIIHHCHHFHHHCEEACCHHHEH)<<CCDGFFDFFFBCB@DFHFHFHHEEFB8@EEEFFEEFFFFFFFFFEFCEFFFCAAC?EF??AC???0*:?C*:::?EE)0>'.42AAECEFE:*0:AAC?D'..8C?:?A)).0001*11::??8A**?A################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:11958:2304 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGGGGAATTTCCTGTGTAGCGGGGAAATGCGTAGATATAGGAAGGAACACCAGGGGCGAAGGCGACCACCGGGACTGAAACTGACACTGAGGGGCGAAAGCGTGGGGGGCAAACG
+
????ABBADDDDDDDDEEEEFFHHHHHIFHHIIIHFEECEFGDECECE5>C:55EEHIHIFGHFGHHHHHFHFFHHFHHHHHFBFEEDEEFFFFEFFFEFEFABEEFFFEEBEFFEFF=::AE*:AEE0?:?CFE8A>'.<?EEE??E?A??CEEF<>'.8AC?ECE)848?0**::AAC???EEE)*0)084'48<'8'882<CA).2<408?*1:??EEE#############################
@MISEQ03:64:000000000-A2H3D:1:1101:19110:2311 1:N:0:TCCACAGGAGT
TACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGGGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTAAAACTGACACTGAGGGGCGAAAGCGGGGGGGGCAAACG
+
????9BB?BBBBBBBBEEECFAFHHHHHHHFHHHHHHCCEFF/CDCCE>DCH5AECFHHHFHHHHHHHHHHHGFHHCFHHHHHEEEDEDEED@EEEEEEEEEEEEEEEEE;;BEEE?EEEEEE?*?CA?EE::?8'.''..?CEE*::/:?A:C?E??82?CCEEEE))4?EEEEA:?*80?AEEC#################################################################
"""

# For reference. Data used to make the 'joined_reads' reference string.
reverse_reads = """@MISEQ03:64:000000000-A2H3D:1:1101:14358:1530 2:N:0:TCCACAGGAGT
ACGGACTACAAGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGGATTTCACATCCAACTTAACAACCCACATACCCGCCTTTTCGCCCAGGTAATCC
+
?????@@BDDBDDDBBCFFCFFHIIIIIIIIFGHHHHEHHHIIIHHHHHFHIIHIGHHIDGGHHHHIIIIICEFHIHHCDEHHHHHHFHHCFHDF?FHHFHHHFFDFFFDEDDD..=DDDE@<BFEEFCFFCECE==CACFE?*0:*CCAA?:*:*:0*A?A80:???A?*00:**0*1*:C??C?A?01*0?);>>'.8::A?###############################################
@MISEQ03:64:000000000-A2H3D:1:1101:14206:1564 2:N:0:TCCACAGGAGT
ACGGACTACAGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGTGCATGAGCGTCAGTGTTAACCCAGGGGGCTGCCTTCGCCATTGGTATTCCTCCACATCTCTACGCATTTCACTGCTACACGTGGAATTCTACCCCCCTCTGCTACACTCTAGCCTTGTAGTTTCAAACGCAGTTCCCAGGTTGAGCCCGGGGCTTTCACATCTGCCTTACAAAACCGCCTGCGCACGCTTTACGCCCCGTAATTC
+
?????@@BDDBDDD?@CFFFFFHHHHHFFHHHHHHHHHHH@FFHEDFFH@FHBGCDHHHBFHHHHHHHEHHHHDCCEFFDFFFEE:=?FF?DFDFDFFF==BEE=DBDDEEEEEB,4??EE@EEE,3,3*3,?:?*0ACCEDD88:***?*0:*0***0*?C?00:AE:?EE:*A8'.?:CAA?A80*0*??AA88;28;C##################################################
@MISEQ03:64:000000000-A2H3D:1:1101:14943:1619 2:N:0:TCCACAGGAGT
ACGGACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCAAGGTTGAGCCCGGGGCTTTCACATCCAACTTACAAAACCACCTACCCGACCTTTACGCCCAGAAATTC
+
?????@@BDDDDDD?AFFFFFFIIHHIHIIHIIIIHIHHHHHHHHHHHHHHHHHHIHHHIIHHIHIIIIII?EFEGHHHHHIIHHDHHFD@FFEHFHFHFHFHFFFFFFFFEEEEFFFDEB<E@EFEEABA=B=CAEFEEFEA?A:*AC::??:**10??:00::*??EC*?:?C*::A*?C*1:8A################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:15764:1660 2:N:0:TCCACAGGAGT
ACGGACTACCGGGGTATCTAATCCTGTTTGCTCCCCTAGCTTTCGCACCTCAGCGTCAGTACCGGGCCAGTGAGCCGCCTTCGCCACTGGTGTTCTTGCGAATATCTACGAATTTCACCTCTACACTCGCAGTTCCACTCACCTCTCCCGGACTCAAGACCGCCAGTACCAAAGGCAGTTCTGGAGTTGAGCTCCAGGTTTTCACCCCTGCTTTAAAAATCCCCCAACGCGGCCTTTCCCCCCAGTGACTC
+
?????@=BB?BBBB<?>ACFFCECFCHCFHH=CGHHDFH=E?ACDEHHCCFFGHHDHH@CBEFHHCHHHF,5@?DF)4<C3D4DD+=,BD5;DE=DBDE=<E<;?E?B;3?;A?;;;EBBC:??EEEEE?;AA*:A??#################################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:15211:1752 2:N:0:TCCACAGGAGT
ACGGACTACCGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGGATTTCCCATACAACTTAACAAACCACTTACGCGGGTTTTCGCCCCACAAATTC
+
?????@@BB-<@BBBBEA?FFFA8ACFCBHHHGHHHHBHHHHFCDDD7CHHFE?CCDDCFGBHHHGGHFGFGFFHFDCDHHC?=CDHFFDCDDDF,EFF5?BFEDBBDB@@EEACCE;,?::@BEEEEACC*??//:AA*8AAAEE?ECC#####################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:15201:1774 2:N:0:TCCACAGGAGT
ACGGACTACCGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTAAGCCCGGGGATTTCCAATCCAACTAAACAAAACACACACCCGCTCTTTACGCCCACCATTTC
+
?????@=BBBBBBB<=CFFFFFFHFCFCEHHDGHHEFEHHHHHHHHHHHFHHGC-AEGFCGHHHFFHFHHBFGDEDDCEDDH+DDDHHF,,7D@DFDFFFBFFEDEDEEDE:@:B?C;,3@<><EEEE*BEEC?E*0AEEEEE*8*:CCE:CA*?*A?:AAA#########################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:15976:1791 2:N:0:TCCACAGGAGT
ACGGACTACAAGGGTTTCTAATCCCGTTTGCTCCCCTAGCTTTCGCACCTCAGCGTCAGAAATGGACCAGAGAGCCGCCTTCGCCACCGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGGTTTTCCAATCCAATTTAACGACCAACCACGCGGCGCTTTAGGCCCGGTAATGC
+
?<???@=BBBBBBB556CFFBCEFFEHHHHHHHHHE=EHHHHHHHHHHHHHHHFHHEDCGFHHHHHHHHH;A?EFHHHHHHHHH+=EHHC)8@+?BFFFDFEEEEE;DDEEEEBCEECEEEBEEEEEEEEEEE:?*/:ACEE)*8*:C:A*0?:A*:C?A:?**:AECE?*?:*:C:????C#####################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:16031:1840 2:N:0:TCCACAGGAGT
ACGGACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTTAACAATACTCTAGCTTGTCAGTTTTAAATGCAGTCCCCAGGTTGAGCCCGGGGCTTTCACATCAAATTTACAAAACCACCCACGCGCGTTCTACGCCGGACAATCC
+
?????==<BBBBBBBBFFFFF?FFF?FFFHHFEFFHHHH:@>CHEDHHHFFFGBCCDDFGGHHHHEHHHHH5AE+C*>==+EDHHDEFCFCDF3@.D=,CFH=,@,4DFFE:=DDDDEB:)1:1;;?B;BE;??,?EE;AEE??**0*/:0??:***:?E*:8?A*:CEE#################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:12964:1852 2:N:0:TCCACAGGAGT
ACGGACTACTCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCAACCACCCTCTACCATACTCTAGCTTGTAAGTTTTGAATGCAGTTCCAAGGTTGAGCCCGGGGCTTTCACACCCAACTTAACAACCCACCTACGCGGCATTTACGCCCACTACTTC
+
?????@=BBBB9?=55AACCC?ECFFFHB>FFHGFHFHHHHHHHHHHHHFHHHGGGHHHGGHHHHHHDDFEGH;EBCEHD+AFE@C+@F=.7D+@CDCFFHHFFFD?DF@E+=:BDDB;D=@BE?BCE*,33,,?,3;;:?).0**::***0/*/0A??:*:****00/**8*0?AE:?AAC**0):??C#############################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17245:1906 2:N:0:TCCACAGGAGT
ACGGACTACAGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGGCTTTCAAATCCAACTTAACAAACAACCAACCGGCGCTTTACGCCCAGTAATTC
+
?????@@BDDDDBDABEFFFFFHHHHHHHHGGHGHHHHHHEH@FEHEEHHHFHHH=EGHHHHDGHHHHFHHGGHHHCCEDEHHHHHHHHHDFHHF=DBDFHFD?BB.BF;@DDDD.=DD*>6)@==ABAACBB5=B,=,88A)???:E*::::??*:**1**8??CCCEE8A:A::AAACAC??A)1:0**1*)48;'42A?EA**1?*1*0::??:ACEF##############################
@MISEQ03:64:000000000-A2H3D:1:1101:18171:1983 2:N:0:TCCACAGGAGT
ACGGACTACCGGGGTTTCTAATCCTGTTCGCTCCCCACGCTTTCGCGCCTGAGCGTCAGGTAAGGCCCAGAGATCCGCCTTCGCCACCGGTGTTCCTCCTGATATCTGCGCATTTCACCGCTACACCAGGAATTCCGATCTCCCCTACATACCTCTAGCCGGGCAGTATCGGATGCAGACCCGAGGTTGAGCCCCGGGATTTCACATCGGCTTACCAAAGCGCCCGGCGCCCCCTTTACGCCCCAGAAACC
+
?????@@BDBDDDD=BCFFFFFIIIIHHFEHHHHIHIHHHEFCGDEHHHEFFEGC>EEHI?EHHGHHHHFH+C=,,?FHDDHFEE@EFE=1;A;EECCE==BEB,BBC=@@<?EE:?E:8>24;:CEAA8?CC*??:0?;*1?AE?CE*10:0*1:CAA*;22;2A#####################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:14225:2009 2:N:0:TCCACAGGAGT
ACGGACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGGCTTTCACATCCAACTTAACAAACAACCACGCCGCGCTTTAGCCCAGGTAATTC
+
?????@@BDDDDDD??CFFFFFHIHHHHIIIIHIIHIHHHHHIHHHHHHFFHHIHHHIHHHHHIIHIHIIIFFFEGHHEDEHHHHDHHHHCFFDFFHHHHHHFFFFFFF@EDEED=DDEED@EBFCEEEFECCCEEEFB<CA888:AEEFEFEA??CCEFF:?:ACCFFE?E:AC?:*:?EFE*:):???::A).;D>D>8:?################################################
@MISEQ03:64:000000000-A2H3D:1:1101:16656:2052 2:N:0:TCCACAGGAGT
ACGGACTACCCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCCCCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCAAGGTTGAGCCCGGGGCTTTCACACCCAACTTAACAAACCACCTACGCGCGCTTTACGCCCAGCATTTC
+
?????@@BDDDDDD<@CFFFFFIHHFHHFHHHHHHHIHHHHHFHCEHHHHIIFHHAFHHHFFHIIHHIIIHGHFEH<DDEDH;DEHHHHHFFH;FFHFHFFFFD?FBFF=BDEDDDFEEAE@BEFFFF<BE=B,=,5?*).;>8A:*:::?E?*::A::?AE8AEFEEEC?A:CE?AEA:EF*008:?EF*:C)8;D228A0:??:*.8A8):*:*1::CE##############################
@MISEQ03:64:000000000-A2H3D:1:1101:18209:2060 2:N:0:TCCACAGGAGT
ACGGACTACTAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGCGTCAATACTTGTCCAGCGAGTCCCCTTCGCCACTGGTGTTCTTCGGAATATCTACGAATTTCACCTCTACACCCGGAATTCCACTCCCCCTTCCAAGATTCCAGCTTAGCAGTTTCAAGGGCAGTTCCGGGGTTGGCACCCGGGATTTCACCCCTGCCTTGCTCAACCCCCCACGGGGCCTTTACCCCCAGCATTCC
+
=9<=9<7<@@@@@@<@A8A>C>8E>-8AE;C99CEEECC>>EECE@CCDE,C@E++5>E-A=E-C@@=@5@C>C<D-5A5CC<CD+<4DE3=<C+4==+<@D++2:9DEE1<9DE########################################################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:13824:2095 2:N:0:TCCACAGGAGT
ACGGACTACCAGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCGTCTCAGCGTCAGTAACGTCCCAGCGAGCTGCCTTCGCAATGGGTGTTCCTCCTGATATCTGCGCATTTCACCGCTACACCAGGAATTCCACTCCCCCCTTCCGTACCCTAGCCCGACAGTACCCCCCGGCTTCCGAGGCTTGACCCCCCGCCTTTCACACCGGACTTACCGGGCCGCCTACCCGGCCTTTCGCCCCACCGTTTC
+
??<??@@BBBBBBBBBCFCFFFHHHHHHBHHGHHHHHCHHEH>GDEHCA:DFGHHEEEEFFHHHHHHDHED7<CHEGHFFDFFHEDHHHDDDE@D??DD;=B,,5B,,56)):?;BEE?*1::):?**8AEAC*?*:?#################################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17839:2106 2:N:0:TCCACAGGAGT
ACGGACTACAAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGAAGTTCCAAGTTTGAGCCCGGGGCTTTCACATCCAATTTAAAAAACAACCAACGCGCGCTTCACCCCAGGTAATAC
+
?????@@BBBBBBB<5ACFFFFHHHHHHHHHHHHHHHHHHHHFFHHHHHHHHHGHHHFHGHHHHFGHHHHH?EFEECCEEEHHHEHHHHHDFHDFDHDHHHHFFDFFFHFEDDDD,5DD@BB<EEEEECBB?B3B;,?,3?E8CE?*?A*/:/:::??AE::**0:AEE##################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17831:2124 2:N:0:TCCACAGGAGT
ACGGACTACAAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTCGAATGCAGTTCCCAGGTTGAGCCCGGGCCTTCAACCTCCACTTTACAAAACAACCAAACGCCGCTTACCGCCACGAAATCC
+
?????@=BBBBBBB5<CFFFFCFHHHHHFHHHHHHHHHHHHHFHEEHHEHFGHGHFGDF?EEFHHHHDGHH=EHEGCCEEEHHHHHH@FFCFH+CFHCFHHHHHFFFHHFE:DDD,=EBDEBE<EEEE?;B?B?E?CEEEEEE8A?CC#######################################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:12555:2129 2:N:0:TCCACAGGAGT
ACGGACTACCGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCCCGTCAGTTTTGAATGCATTTCCCAGGTTGAGCCCGGGGCTTTCAAACCCAACTAAACAAACAACCAACGCGCGTTTTCCGCCACGTAATTC
+
?????@@BDDDDDDDDEFFFFFHHHFHIIIHHHIIHHHHHEHFHEHHHHIFGHHIHIHFHI=FHFIHIIIHDFHHHHEHEDHHHHHHHHHHHH@FBFFFFFFFEFDEFE6@:@ACBEFFEEB@BCB=A<<A?C:A::C8AEE)0?A*?A:*:**10?1**1.4A88?EE?#################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:13627:2170 2:N:0:TCCACAGGAGT
ACGGACTACTAGGGTATCTAATCCCGTTTGCTACCCTAGCTTTCGCACCTCAGTGTCAGATTGGGTCCAGTGAGCCGCTTTCGCCACCGGTGTTCCTCCAGATATCTACGCATTTCACCGCTCCACCTGGAGTTCCGCTCACCCCTACCCACCTCGAGCCGGCCAGTATCCCGCGCAATTCCCCGTTTGAGCCGGGGGTTTTCCCAAGGGTCTTAACAGCCCACCTACGTCCCCTTTATGCCAGGTAATTC
+
?????@=BD?BBDD?58ACFFCHHHHHHHHFGHHHEEFHHHHHCDEEEECFDGFHDGHHHHFFFHEHHHHHHHFFFHAEHFFEHHHHFHH<DE:C--@F-CCF=,,4BDFE:@E,BBED@:)2>=C?;BC=?C,==*3.84?EC?88A8A:A?*8?###############################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:11781:2223 2:N:0:TCCACAGGAGT
ACGGACTACTGGGGTATCTAATCCTGTTCGCTCCCCATGCTTTCGCTCCTCAGCGTCAGTTAATGCCCAGAGACCTGCCTTCGCCATCGGTGTTCCTCCTGATATCTGCGCATTTCACCGCTACACCAGGAATTCCAATCTCCCCTACATCACTCTAGTCTGCCCGTACCCACCGCAGATCCGGAGTTGAGCCCCGGACTTTCACGGCAGACGCGCAAACCGCCCACAGAGCTCTTTCCCCCCAAAAATCC
+
?????@@BDDDDDDBDCFFFFFHHHHHHHHHHHEHEHDFHHHHHHEDEEHHHFHHHHEHHHHHFGHHHHDD;EFFHCFECAGFEEE+@E@3?E:DFFFHBDHC?BBDFFE8;=DD,,=DEE<@),==?B*44=,=,4**0*0:CA*A?::*0::?0AC:CE*?:*8.4AE?*8?)'82;*?0?####################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:17996:2254 2:N:0:TCCACAGGAGT
ACGGACTACCAGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTGGAGCCCGGGGCTTTACAATCAAACTTAACAAACCCCCTACGCGCGCTTTAGCCCCACAATTTC
+
?9???@@BBDDDDD<BEFFFFFHHHHCEFHGGHHHHHAEHEHFHEEHHCBFHFFHEGHFCGGFCGHHFHHFGHHHHDHHECEDEHDDHH?@?FDFDHHHFHFFDDHDFFFEEEEEDEDB=>BBEE=BECEEEE,B=?CBACE)*0**?C?*:*0:*:?:A:??**::?E**80::::??:CAC:C8C################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:13712:2276 2:N:0:TCCACAGGAGT
ACGGACTACTGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGAAGTTCCAAGGTTAAGCCGGGGGCTTTCACATCCAACTAACCAAACCACTAACCGCCGCTTTAGGCCCAGCAATTC
+
?????@@BDDDDDDBDEFFFFFHFHHHIFFHIHIIHHEHHHHFHHHEHCFHHDFGGHIIHIHFGHHHGFHF-AEEGHHHHEGFHHHDFHB?FHEHCF?FDDFH??BFFFF=DDEEEFFEDE8>:BECCAEEAECE,ABAACEA*0**?*01?001*::*A**0:*::CE1*8:?**11:*CE#####################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:15819:2301 2:N:0:TCCACAGGAGT
ACGGACTACCGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGAAGTTCCCAGGTTGACCCCGGGGCTTTCACAGCAGCCTTAACTCACCCCCTGCCAACCCTTTACTCCCCGAAATTC
+
?9???@=BDDDDDD<@CFFFFFHHHBFFHHHFGHHHHHHHEHFGEHHHHEDGD?FCGHHFHHHHHHHHBDF?FHHHFHH@DHHHH+DHHDCFHDFDFBFFDFFEDFEEDEEDEEC=CCEEEBEEFCEFEEFEE?:CEE*8CA800*:E*:AA::***1??*:**1::CF##################################################################################
@MISEQ03:64:000000000-A2H3D:1:1101:11958:2304 2:N:0:TCCACAGGAGT
ACGGACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGCCTTTAAAATCCAACTTAAAAAACAACTAACCGCCGCTTTCCGCCCAGAAATCC
+
?????@@BDBDDD@@BFFFFFFHFHHHHFHHHHHHHHHHHHHHHHHHHHHFHHHHHFGH=CGEH=FHFGHFEFHEHHCEEEHHHDCEFHH.?FDFFHHHHBFHFFHFFFFEEEEEEEEEEB@EEEECE;CC?BCEEEC;;CEA*8:AE**00::C0A::?:*0*AEE?E?*A**:C?*:?:?**0)00::**?8>'8';ACA*0*0C?:?******::??E?CE###########################
@MISEQ03:64:000000000-A2H3D:1:1101:19110:2311 2:N:0:TCCACAGGAGT
ACGGACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGGGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTAAGCCCGGGGCTTTCCAACCCAACTTAACAAACCACCTACCCCCGCTTTACGCCCAGAAATTC
+
?????@=BDDDDDD<5CFFFFFHHHHHF>CGFHHHHHHHEEHFHEHHHHHHHGAFGHHHHH-5AF5AFHHD+5*5CCCDDHFFHEEHFHHBFF:BDD4?=.=DEFFEBEDEBEEECEFFEE<::CEAACE?:A*1:*C88?AE.?:*::**1:**11*01***1?C??:?0?:C:C**1*1::*:8A?10*1?##########################################################
"""

if __name__ == "__main__":
    main()
