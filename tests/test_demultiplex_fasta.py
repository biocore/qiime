#!/usr/bin/env python
# file test_demultiplex_fasta.py

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "William Walters"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"

from os import close
from os.path import exists, join
from shutil import rmtree
from tempfile import mkdtemp, mkstemp

from skbio.util import create_dir
from unittest import TestCase, main
from skbio.util import remove_files

from qiime.parse import parse_qual_score
from qiime.demultiplex_fasta import (
    demultiplex_sequences, assign_seqs, check_map, get_ids_bcs_added_field,
    get_bc_lens, process_log_data, process_bc_freqs,
    get_added_demultiplex_field, get_exact_bc_matches, attempt_bc_correction,
    get_curr_bc_added_field, get_demultiplex_data, write_qual_line,
    write_fasta_line, get_label_line, initialize_log_data,
    get_output_ids, assign_seqs, process_files_and_demultiplex_sequences
)


class FakeOutFile(object):

    def __init__(self, name="test_file"):
        self.data = ""
        self.name = name

    def write(self, s):
        self.data += s


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """ Set up test data, output dir for unit tests """

        self.valid_mapping_data_golay = valid_mapping_data_golay
        self.valid_mapping_data_golay_upper = valid_mapping_data_golay_upper
        self.valid_mapping_data_var_len_bcs = valid_mapping_data_var_len_bcs
        self.valid_mapping_data_no_bcs_added_demultiplex =\
            valid_mapping_data_no_bcs_added_demultiplex
        self.valid_mapping_data_bcs_and_added_demultiplex =\
            valid_mapping_data_bcs_and_added_demultiplex
        self.invalid_mapping_data_golay_header =\
            invalid_mapping_data_golay_header
        self.invalid_mapping_data_golay_dna =\
            invalid_mapping_data_golay_dna
        self.invalid_mapping_data_golay_dup_bcs =\
            invalid_mapping_data_golay_dup_bcs
        self.invalid_mapping_data_golay_dup_sids =\
            invalid_mapping_data_golay_dup_sids
        self.invalid_mapping_data_no_bcs_added_demultiplex =\
            invalid_mapping_data_no_bcs_added_demultiplex
        self.invalid_mapping_data_bcs_added_demultiplex =\
            invalid_mapping_data_bcs_added_demultiplex
        self.invalid_mapping_data_golay_missing_bc =\
            invalid_mapping_data_golay_missing_bc
        self.expected_formatted_log_data_no_qual =\
            expected_formatted_log_data_no_qual
        self.expected_formatted_log_data_with_qual =\
            expected_formatted_log_data_with_qual
        self.valid_mapping_data_no_bcs_no_added_demultiplex =\
            valid_mapping_data_no_bcs_no_added_demultiplex
        self.valid_fasta_file_no_errors = valid_fasta_file_no_errors
        self.valid_qual_file_no_errors = valid_qual_file_no_errors
        self.valid_fasta_file_with_bc_errors = valid_fasta_file_with_bc_errors
        self.sample_correct_mapping_data = sample_correct_mapping_data
        self.sample_fasta_file = sample_fasta_file
        self.sample_qual_file = sample_qual_file

        self.output_dir = mkdtemp()
        self.output_dir += '/'
        create_dir(self.output_dir)

        fd, self.correct_mapping_fp = mkstemp(
            prefix='correct_mapping_',
            suffix='.txt')
        close(fd)
        map_file = open(self.correct_mapping_fp, 'w')
        map_file.write(self.sample_correct_mapping_data)
        map_file.close()

        fd, self.sample_fasta_fp = mkstemp(
            prefix='sample_fasta_',
            suffix='.fna')
        close(fd)
        sample_fasta = open(self.sample_fasta_fp, 'w')
        sample_fasta.write(self.sample_fasta_file)
        sample_fasta.close()

        fd, self.sample_qual_fp = mkstemp(
            prefix='sample_qual_',
            suffix='.qual')
        close(fd)
        sample_qual = open(self.sample_qual_fp, 'w')
        sample_qual.write(self.sample_qual_file)
        sample_qual.close()

        self._files_to_remove =\
            [self.correct_mapping_fp, self.sample_fasta_fp,
             self.sample_qual_fp]

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)

    def test_process_files_and_demultiplex_sequences(self):
        """ Overall IO/functionality test """

        process_files_and_demultiplex_sequences(self.correct_mapping_fp,
                                                [self.sample_fasta_fp], [self.sample_qual_fp],
                                                output_dir=self.output_dir, write_unassigned_reads=True,
                                                save_barcode_frequencies=True)

        log_file = open(join(self.output_dir, "demultiplex_fasta.log"), "U")
        log_lines = [line for line in log_file]
        log_file.close()

        bc_freqs_f = open(join(self.output_dir, "barcode_freqs.txt"), "U")
        bc_freqs = set([line for line in bc_freqs_f])
        bc_freqs_f.close()

        fasta_seqs_f = open(join(self.output_dir, "demultiplexed_seqs.fna"),
                            "U")
        fasta_lines = [line for line in fasta_seqs_f]
        fasta_seqs_f.close()

        qual_seqs_f = open(join(self.output_dir, "demultiplexed_seqs.qual"),
                           "U")
        qual_lines = [line for line in qual_seqs_f]
        qual_seqs_f.close()

        unassigned_fna_f = open(join(self.output_dir, "unassigned_seqs.fna"),
                                "U")
        unassigned_fna_lines = [line for line in unassigned_fna_f]
        unassigned_fna_f.close()

        unassigned_qual_f = open(join(self.output_dir, "unassigned_seqs.qual"),
                                 "U")
        unassigned_qual_line = [line for line in unassigned_qual_f]
        unassigned_qual_f.close()

        expected_log_lines = [
            'Total sequences in input files:\t4\n', 'Retain barcode:\tFalse\n',
            'Barcode type:\tgolay_12\n',
            'Max barcode error/mismatches allowed:\t0.5\n',
            'Starting sequence identifier:\t1\n',
            'Write unassigned reads:\tTrue\n',
            'Disable barcode correction:\tFalse\n',
            'Added demultiplex field:\tNone\n',
            'Save barcode frequencies:\tFalse\n', '\n',
            'Barcodes corrected/not corrected:\t0/1\n',
            'Number of samples in mapping file:\t3\n',
            'Sample count min/max/mean:\t1 / 1 / 1.00\n',
            'Sample\tSequence Count\tBarcode/Added Demultiplex\n',
            'PC.354\t1\tAGCACGAGCCTA\n', 'PC.356\t1\tACAGACCACTCA\n',
            'PC.355\t1\tAACTCGTCGATG\n', 'Seqs written\t3\n',
            'Percent of input seqs written\t0.75']

        # Use set() to avoid order issues with unit tests.
        expected_bc_freqs = set(['Barcode frequencies\n',
                                 'ACCCCCCACTCA\t1\n', 'AACTCGTCGATG\t1\n',
                                 'AGCACGAGCCTA\t1\n', 'ACAGACCACTCA\t1'])

        expected_fasta_lines = [
            '>PC.354_1 ABCD0001 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0\n',
            'CAGGACGAGACGAGGTT\n', '>PC.355_2 EFGH0002 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n',
            'CCAGATTACGAGATTA\n', '>PC.356_3 IJKL0003 orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0\n',
            'GACCGATTACGATAACG\n']

        expected_qual_lines = [
            '>PC.354_1 ABCD0001 orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0\n',
            '30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n',
            '>PC.355_2 EFGH0002 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n',
            '12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n',
            '>PC.356_3 IJKL0003 orig_bc=ACAGACCACTCA new_bc=ACAGACCACTCA bc_diffs=0\n',
            '10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n']

        expected_unassigned_fna_lines = [
            '>Unassigned_4 nomatchingBC orig_bc=ACCCCCCACTCA new_bc=GTCCCCCACTGA bc_diffs=3\n',
            'ACCCCCCACTCAGACCGATTACGATAACG\n']

        expected_unassigned_qual_line = [
            '>Unassigned_4 nomatchingBC orig_bc=ACCCCCCACTCA new_bc=GTCCCCCACTGA bc_diffs=3\n',
            '30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n']

        # Skips first few lines that point to tmp filepaths
        self.assertEqual(log_lines[5:], expected_log_lines)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(fasta_lines, expected_fasta_lines)
        self.assertEqual(qual_lines, expected_qual_lines)
        self.assertEqual(unassigned_fna_lines, expected_unassigned_fna_lines)
        self.assertEqual(unassigned_qual_line, expected_unassigned_qual_line)

    def test_check_map(self):
        """ Properly returns data with valid input """

        header, mapping_data = check_map(self.valid_mapping_data_golay)

        expected_header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        expected_mapping_data =\
            [['s1', 'AACTCGTCGATG', 'ATTCGATART', 's1_description'],
             ['s2', 'agcAGCACTTGT', 'ATTCGATART', 's2_description'],
                ['s3', 'ACCGCAGAGTCA', 'YATGCTGCCTCCCGTAGGAGT', 's3_description']]

        self.assertEquals(header, expected_header)
        self.assertEquals(mapping_data, expected_mapping_data)

    def test_check_map_errors(self):
        """ Properly raises errors with invalid input """

        # Bad header
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_golay_header)
        # non DNA characters
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_golay_dna)
        # Duplicate barcodes
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_golay_dup_bcs)
        # Duplicate SampleIDs
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_golay_dup_sids)
        # More than one SampleID, no barcodes or added demultiplex specified
        self.assertRaises(ValueError, check_map,
                          self.valid_mapping_data_no_bcs_added_demultiplex, barcode_type=0)
        # No barcodes, added_demultiplex has duplicates
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_no_bcs_added_demultiplex, barcode_type=0,
                          added_demultiplex_field="Added_Demultiplex")
        # Barcodes plus added demultiplex results in duplications
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_bcs_added_demultiplex,
                          added_demultiplex_field="Added_Demultiplex")
        # Missing a barcode
        self.assertRaises(ValueError, check_map,
                          self.invalid_mapping_data_golay_missing_bc,
                          barcode_type="variable_length")

    def test_check_map_var_len_bcs(self):
        """ Properly handles variable length barcodes """

        header, mapping_data = check_map(self.valid_mapping_data_var_len_bcs,
                                         barcode_type="variable_length")

        expected_header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        expected_mapping_data =\
            [['s1', 'AACTCGTCGATG', 'ATTCGATART', 's1_description'],
             ['s2', 'CACTTGT', 'ATTCGATART', 's2_description'],
                ['s3', 'ACCGCAGAGTCA', 'YATGCTGCCTCCCGTAGGAGT', 's3_description']]

        self.assertEquals(header, expected_header)
        self.assertEquals(mapping_data, expected_mapping_data)

    def test_check_map_var_len_not_specified(self):
        """ Raises error if var len bcs detected but not specified """

        self.assertRaises(ValueError, check_map,
                          self.valid_mapping_data_var_len_bcs)

    def test_check_map_added_demultiplex(self):
        """ Properly handles no barcodes plus added_demultiplex field """

        header, mapping_data = check_map(
            self.valid_mapping_data_no_bcs_added_demultiplex,
            barcode_type=0, added_demultiplex_field="Added_Demultiplex")

        expected_header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
             'Added_Demultiplex', 'Description']
        expected_mapping_data =\
            [['s1', '', 'ATTCGATART', '1', 's1_description'],
             ['s2', '', 'ATTCGATART', '2', 's2_description'],
                ['s3', '', 'YATGCTGCCTCCCGTAGGAGT', '3', 's3_description']]

        self.assertEquals(header, expected_header)
        self.assertEquals(mapping_data, expected_mapping_data)

    def test_check_map_single_sample(self):
        """ Properly handles single sample with no demultiplex data"""

        header, mapping_data = check_map(
            valid_mapping_data_no_bcs_no_added_demultiplex,
            barcode_type=0)

        expected_header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        expected_mapping_data =\
            [['s1', '', '', 's1_description']]

        self.assertEquals(header, expected_header)
        self.assertEquals(mapping_data, expected_mapping_data)

    def test_check_map_bcs_and_added_demultiplex(self):
        """ Properly handles barcodes plus added_demultiplex field """

        # Combinations of bc plus added demultiplex field need to be unique
        # but can be duplicated if both options are being used.
        header, mapping_data = check_map(
            self.valid_mapping_data_bcs_and_added_demultiplex,
            barcode_type=4, added_demultiplex_field="Added_Demultiplex")

        expected_header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
             'Added_Demultiplex', 'Description']
        expected_mapping_data =\
            [['s1', 'AAAA', 'ATTCGATART', '1', 's1_description'],
             ['s2', 'AAAA', 'ATTCGATART', '2', 's2_description'],
                ['s3', 'CCCC', 'YATGCTGCCTCCCGTAGGAGT', '1', 's3_description']]

        self.assertEquals(header, expected_header)
        self.assertEquals(mapping_data, expected_mapping_data)

    def test_get_ids_bcs_added_field(self):
        """ Properly returns list of lists of SampleID/bcs/demultiplex data """

        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'Added_Demultiplex', 'Description']
        mapping_data = [['s1', 'AAAA', 'ATTCGATART', '1', 's1_description'],
                        ['s2', 'ttTT', 'ATTCGATART', '2', 's2_description'],
                        ['s3', 'CCCC', 'YATGCTGCCTCCCGTAGGAGT', '3', 's3_description']]

        actual_ids_bcs_added_field = get_ids_bcs_added_field(header,
                                                             mapping_data)

        expected_ids_bcs_added_field =\
            {('CCCC', ''): 's3', ('AAAA', ''): 's1', ('TTTT', ''): 's2'}

        self.assertEqual(actual_ids_bcs_added_field,
                         expected_ids_bcs_added_field)

    def test_get_ids_bcs_added_demultiplex_field(self):
        """ Properly returns list of lists of SampleID/bcs/demultiplex data """

        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'Added_Demultiplex', 'Description']
        mapping_data = [['s1', 'AAAA', 'ATTCGATART', '1', 's1_description'],
                        ['s2', 'TTTT', 'ATTCGATART', '2', 's2_description'],
                        ['s3', 'CCCC', 'YATGCTGCCTCCCGTAGGAGT', '3', 's3_description']]

        actual_ids_bcs_added_field = get_ids_bcs_added_field(header,
                                                             mapping_data, added_demultiplex_field='Added_Demultiplex')

        expected_ids_bcs_added_field =\
            {('AAAA', '1'): 's1', ('CCCC', '3'): 's3', ('TTTT', '2'): 's2'}

        self.assertEqual(actual_ids_bcs_added_field,
                         expected_ids_bcs_added_field)

    def test_get_ids_bcs_no_bcs(self):
        """ Properly returns list of lists of SampleID/bcs/demultiplex data """

        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'Added_Demultiplex', 'Description']
        mapping_data = [['s1', 'AAAA', 'ATTCGATART', '1', 's1_description'],
                        ['s2', 'TTTT', 'ATTCGATART', '2', 's2_description'],
                        ['s3', 'CCCC', 'YATGCTGCCTCCCGTAGGAGT', '3', 's3_description']]

        actual_ids_bcs_added_field = get_ids_bcs_added_field(header,
                                                             mapping_data, barcode_type=0,
                                                             added_demultiplex_field='Added_Demultiplex')

        expected_ids_bcs_added_field =\
            {('', '2'): 's2', ('', '1'): 's1', ('', '3'): 's3'}

        self.assertEqual(actual_ids_bcs_added_field,
                         expected_ids_bcs_added_field)

    def test_get_bc_lens(self):
        """ Properly returns lengths of barcodes """

        sample_data = {('CCCC', ''): 's3', ('AAAA', ''): 's1',
                       ('TTTT', ''): 's2'}

        expected_lens = [4]

        actual_lens = get_bc_lens(sample_data)

        self.assertEqual(actual_lens, expected_lens)

        # Test with multiple lengths

        sample_data = {('CCCC', ''): 's3', ('', ''): 's1',
                       ('TTTTT', ''): 's2'}

        expected_lens = [5, 4, 0]

        actual_lens = get_bc_lens(sample_data)

        self.assertEqual(actual_lens, expected_lens)

    def test_process_log_data_no_qual(self):
        """ Properly processes and formats log data with no qual files """

        log_data = {'AAAA': 10, 'TTTT': 20, 'CCCC': 0}
        seq_counts = 32
        mapping_file = 'test_mapping.txt'
        fasta_files = [FakeOutFile(name="fasta1.fna"),
                       FakeOutFile(name="fasta2.fna")]
        qual_files = []
        corrected_bc_count = [7, 15]

        actual_log_data = process_log_data(log_data, seq_counts, mapping_file,
                                           fasta_files, qual_files, corrected_bc_count)

        expected_log_data = self.expected_formatted_log_data_no_qual

        self.assertEqual(actual_log_data, expected_log_data)

    def test_process_log_data_with_qual(self):
        """ Properly processes and formats log data with qual files """

        log_data = {'AAAA': 10, 'TTTT': 20, 'CCCC': 0}
        seq_counts = 32
        mapping_file = 'test_mapping.txt'
        fasta_files = [FakeOutFile(name="fasta1.fna")]
        qual_files = [FakeOutFile(name="qual1.qual")]
        corrected_bc_count = [7, 15]

        actual_log_data = process_log_data(log_data, seq_counts, mapping_file,
                                           fasta_files, qual_files, corrected_bc_count, keep_barcode=True,
                                           barcode_type=8, max_bc_errors=0, start_index=1000,
                                           write_unassigned_reads=True, disable_bc_correction=True,
                                           added_demultiplex_field='run_prefix', save_barcode_frequencies=True)

        expected_log_data = self.expected_formatted_log_data_with_qual

        self.assertEqual(actual_log_data, expected_log_data)

    def test_process_bc_freqs(self):
        """ Properly sorts barcode frequency from largest to smallest """

        bc_freqs = {'AAAA': 20, 'TTTT': 47, 'CCCC': 5, 'GGGG': 0}

        actual_freqs = process_bc_freqs(bc_freqs)

        expected_freqs = ['TTTT\t47', 'AAAA\t20', 'CCCC\t5', 'GGGG\t0']

        self.assertEqual(actual_freqs, expected_freqs)

    def test_get_added_demultiplex_field_run_prefix(self):
        """ Properly returns added demultiplex data using run_prefix """

        ids_bcs_added_field = {("AAAA", "ABC"): "S1", ("TTTT", "239A"): "S2"}
        fasta_label = "ABC1234 region=1 length=229"
        added_demultiplex_field = "run_prefix"
        actual_result = get_added_demultiplex_field(ids_bcs_added_field,
                                                    fasta_label, added_demultiplex_field)

        expected_result = "ABC"
        self.assertEqual(actual_result, expected_result)

    def test_get_added_demultiplex_field_run_prefix_no_matches(self):
        """ Properly returns None when no match with run_prefix """

        ids_bcs_added_field = {("AAAA", "ABC"): "S1", ("TTTT", "239A"): "S2"}
        fasta_label = "xxABC1234 region=1 length=229"
        added_demultiplex_field = "run_prefix"
        actual_result = get_added_demultiplex_field(ids_bcs_added_field,
                                                    fasta_label, added_demultiplex_field)

        expected_result = None
        self.assertEqual(actual_result, expected_result)

    def test_get_added_demultiplex_field(self):
        """ Properly returns added demultiplex data using field in line """

        ids_bcs_added_field = {("AAAA", "2"): "S1", ("TTTT", "1"): "S2"}
        fasta_label = "ABC1234 region=1 length=229"
        added_demultiplex_field = "region"
        actual_result = get_added_demultiplex_field(ids_bcs_added_field,
                                                    fasta_label, added_demultiplex_field)

        expected_result = "1"
        self.assertEqual(actual_result, expected_result)

    def test_get_added_demultiplex_field_no_match(self):
        """ Properly returns None with no matches to field in line """

        ids_bcs_added_field = {("AAAA", "2"): "S1", ("TTTT", "3"): "S2"}
        fasta_label = "ABC1234 region=1 length=229"
        added_demultiplex_field = "region"
        actual_result = get_added_demultiplex_field(ids_bcs_added_field,
                                                    fasta_label, added_demultiplex_field)

        expected_result = None
        self.assertEqual(actual_result, expected_result)

    def test_get_exact_bc_matches_hit(self):
        """ Returns exact matches to barcodes """

        curr_bc = "AAAA"
        all_bcs = ["AAAAA", "TTTT", "CC", "AAAA"]
        actual_bc = get_exact_bc_matches(curr_bc, all_bcs)

        expected_bc = "AAAA"
        self.assertEqual(actual_bc, expected_bc)

    def test_get_exact_bc_matches_no_hit(self):
        """ Returns None if no barcodes match """

        curr_bc = "AAAA"
        all_bcs = ["AAAAA", "TTTT", "CC"]
        actual_bc = get_exact_bc_matches(curr_bc, all_bcs)

        expected_bc = None
        self.assertEqual(actual_bc, expected_bc)

    def test_attempt_bc_correction_exact_match(self):
        """ Corrects barcode with exact match to given barcodes """

        curr_bc = "AGCAGCACTTGT"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12"
        actual_bc, actual_errs = attempt_bc_correction(curr_bc,
                                                       all_bcs, barcode_type)

        expected_bc = "AGCAGCACTTGT"
        expected_errs = 0
        self.assertEqual(actual_bc, expected_bc)
        self.assertEqual(actual_errs, expected_errs)

    def test_attempt_bc_correction_golay12(self):
        """ Corrects barcode using golay correction """

        curr_bc = "GGCAGCACTTGT"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12"
        actual_bc, actual_errs = attempt_bc_correction(curr_bc,
                                                       all_bcs, barcode_type)

        expected_bc = "AGCAGCACTTGT"
        expected_errs = 1
        self.assertEqual(actual_bc, expected_bc)
        self.assertEqual(actual_errs, expected_errs)

    def test_attempt_bc_correction_golay12_exceeds_range(self):
        """ Returns None, count of errors when exceeds correction range """

        # Has 4 changes from a valid golay code
        curr_bc = "TCGTGCACTTGT"
        all_bcs = ["AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12"
        actual_bc, actual_errs = attempt_bc_correction(curr_bc,
                                                       all_bcs, barcode_type)

        expected_bc = None
        expected_errs = 4
        self.assertEqual(actual_bc, expected_bc)
        self.assertEqual(actual_errs, expected_errs)

    def test_attempt_bc_correction_hamming8(self):
        """ Corrects barcode using hamming correction """

        curr_bc = "AGCAGCAC"
        all_bcs = ["AGCAGAAC", "TGCAGTAC", "ACAGAGTC"]
        barcode_type = "hamming_8"
        actual_bc, actual_errs = attempt_bc_correction(curr_bc,
                                                       all_bcs, barcode_type)

        expected_bc = "AGCAGAAC"
        expected_errs = 0.5
        self.assertEqual(actual_bc, expected_bc)
        self.assertEqual(actual_errs, expected_errs)

    def test_attempt_bc_correction_generic(self):
        """ Corrects barcode using generic correction """

        curr_bc = "GGCAGCACTA"
        all_bcs = ["AACTCGTCGA", "AGCAGCACTT", "ACAGAGTCGG"]
        barcode_type = 10
        actual_bc, actual_errs = attempt_bc_correction(curr_bc,
                                                       all_bcs, barcode_type)

        expected_bc = "AGCAGCACTT"
        expected_errs = 2
        self.assertEqual(actual_bc, expected_bc)
        self.assertEqual(actual_errs, expected_errs)

    def test_attempt_bc_correction_no_barcode(self):
        """ Returns empty string, no errors for zero length barcode """

        curr_bc = ""
        all_bcs = [""]
        barcode_type = 0
        actual_bc, actual_errs = attempt_bc_correction(curr_bc,
                                                       all_bcs, barcode_type)

        expected_bc = ""
        expected_errs = 0
        self.assertEqual(actual_bc, expected_bc)
        self.assertEqual(actual_errs, expected_errs)

    def test_get_curr_bc_added_field_barcode_only(self):
        """ Gets corrected barcode, no added demultiplex field """

        curr_bc = "GGCAGCACTTGT"
        ids_bcs_added_field = {("AACTCGTCGATG", ""): "s1",
                               ("AGCAGCACTTGT", ""): "s2", ("ACAGAGTCGGCT", ""): "s3"}
        fasta_label = "123ABC region=1 length=255"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12",
        disable_bc_correction = False
        added_demultiplex_field = None
        corrected_bc, num_errors, added_field =\
            get_curr_bc_added_field(curr_bc, ids_bcs_added_field, fasta_label,
                                    all_bcs, barcode_type, disable_bc_correction, added_demultiplex_field)

        expected_corrected_bc = "AGCAGCACTTGT"
        expected_num_errors = 1
        expected_added_field = None

        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_curr_bc_added_field_barcode_only_disabled_correction(self):
        """ Returns None when no exact match and correction disabled """

        curr_bc = "GGCAGCACTTGT"
        ids_bcs_added_field = {("AACTCGTCGATG", ""): "s1",
                               ("AGCAGCACTTGT", ""): "s2", ("ACAGAGTCGGCT", ""): "s3"}
        fasta_label = "123ABC region=1 length=255"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12",
        disable_bc_correction = True
        added_demultiplex_field = None
        corrected_bc, num_errors, added_field =\
            get_curr_bc_added_field(curr_bc, ids_bcs_added_field, fasta_label,
                                    all_bcs, barcode_type, disable_bc_correction, added_demultiplex_field)

        expected_corrected_bc = None
        expected_num_errors = 0
        expected_added_field = None

        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_curr_bc_added_field_barcode_only_no_errs_no_correction(self):
        """ Returns barcode with exact match, no error correction """

        curr_bc = "AGCAGCACTTGT"
        ids_bcs_added_field = {("AACTCGTCGATG", ""): "s1",
                               ("AGCAGCACTTGT", ""): "s2", ("ACAGAGTCGGCT", ""): "s3"}
        fasta_label = "123ABC region=1 length=255"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12",
        disable_bc_correction = True
        added_demultiplex_field = None
        corrected_bc, num_errors, added_field =\
            get_curr_bc_added_field(curr_bc, ids_bcs_added_field, fasta_label,
                                    all_bcs, barcode_type, disable_bc_correction, added_demultiplex_field)

        expected_corrected_bc = "AGCAGCACTTGT"
        expected_num_errors = 0
        expected_added_field = None

        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_curr_bc_added_field_bc_and_added_field(self):
        """ Gets corrected barcode and added demultiplex field data """

        curr_bc = "GGCAGCACTTGT"
        ids_bcs_added_field = {("AACTCGTCGATG", "123"): "s1",
                               ("AGCAGCACTTGT", "123"): "s2", ("ACAGAGTCGGCT", "ABCD"): "s3"}
        fasta_label = "123ABC region=1 length=255"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12",
        disable_bc_correction = False
        added_demultiplex_field = 'run_prefix'
        corrected_bc, num_errors, added_field =\
            get_curr_bc_added_field(curr_bc, ids_bcs_added_field, fasta_label,
                                    all_bcs, barcode_type, disable_bc_correction, added_demultiplex_field)

        expected_corrected_bc = "AGCAGCACTTGT"
        expected_num_errors = 1
        expected_added_field = "123"

        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_curr_bc_added_field_bc_and_added_field_no_hit(self):
        """ Gets corrected barcode and no matches to added demultiplex field"""

        curr_bc = ""
        ids_bcs_added_field = {("", "1"): "s1",
                               ("", "2"): "s2", ("", "3"): "s3"}
        fasta_label = "123ABC region=1 length=255"
        all_bcs = [""]
        barcode_type = 0,
        disable_bc_correction = False
        added_demultiplex_field = 'region'
        corrected_bc, num_errors, added_field =\
            get_curr_bc_added_field(curr_bc, ids_bcs_added_field, fasta_label,
                                    all_bcs, barcode_type, disable_bc_correction, added_demultiplex_field)

        expected_corrected_bc = ""
        expected_num_errors = 0
        expected_added_field = "1"

        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_curr_bc_added_field_added_field_only(self):
        """ Gets data with no barcode, added field only """

        curr_bc = "GGCAGCACTTGT"
        ids_bcs_added_field = {("AACTCGTCGATG", "123"): "s1",
                               ("AGCAGCACTTGT", "123"): "s2", ("ACAGAGTCGGCT", "ABCD"): "s3"}
        fasta_label = "123ABC region=1 length=255"
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12",
        disable_bc_correction = False
        added_demultiplex_field = 'run_prefix'
        corrected_bc, num_errors, added_field =\
            get_curr_bc_added_field(curr_bc, ids_bcs_added_field, fasta_label,
                                    all_bcs, barcode_type, disable_bc_correction, added_demultiplex_field)

        expected_corrected_bc = "AGCAGCACTTGT"
        expected_num_errors = 1
        expected_added_field = "123"

        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_demultiplex_data_fixed_length_bcs(self):
        """ Gets bc/added data, handles fixed length barcodes """

        ids_bcs_added_field = {("AACTCGTCGATG", ""): "s1",
                               ("AGCAGCACTTGT", ""): "s2", ("ACAGAGTCGGCT", ""): "s3"}
        fasta_label = "123ABC region=1 length=255"
        fasta_seq = "TGCAGCACTTGTACAGATTAGACCGAGATTACGAGATTACCAAGAT"
        bc_lens = [12]
        all_bcs = ["AACTCGTCGATG", "AGCAGCACTTGT", "ACAGAGTCGGCT"]
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        disable_bc_correction = False
        added_demultiplex_field = None
        curr_bc, corrected_bc, num_errors, added_field =\
            get_demultiplex_data(ids_bcs_added_field, fasta_label, fasta_seq,
                                 bc_lens, all_bcs, barcode_type, max_bc_errors, disable_bc_correction,
                                 added_demultiplex_field)

        expected_curr_bc = "TGCAGCACTTGT"
        expected_corrected_bc = "AGCAGCACTTGT"
        expected_num_errors = 1
        expected_added_field = None

        self.assertEqual(curr_bc, expected_curr_bc)
        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_demultiplex_data_variable_length_bcs(self):
        """ Gets bc/added data, handles variable length barcodes """

        # Should return on first hit, from largest to shortest barcode
        ids_bcs_added_field = {("AACTCGT", ""): "s1",
                               ("AGCAGCACTTGT", ""): "s2", ("AGCAG", ""): "s3"}
        fasta_label = "123ABC region=1 length=255"
        fasta_seq = "AGCAGCACTTGTACAGATTAGACCGAGATTACGAGATTACCAAGAT"
        bc_lens = [12, 7, 5]
        all_bcs = ["AACTCGT", "AGCAGCACTTGT", "AGCAG"]
        barcode_type = "variable_length"
        max_bc_errors = 1.5
        disable_bc_correction = True
        added_demultiplex_field = None
        curr_bc, corrected_bc, num_errors, added_field =\
            get_demultiplex_data(ids_bcs_added_field, fasta_label, fasta_seq,
                                 bc_lens, all_bcs, barcode_type, max_bc_errors, disable_bc_correction,
                                 added_demultiplex_field)

        expected_curr_bc = "AGCAGCACTTGT"
        expected_corrected_bc = "AGCAGCACTTGT"
        expected_num_errors = 0
        expected_added_field = None

        self.assertEqual(curr_bc, expected_curr_bc)
        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_get_demultiplex_data_variable_length_bcs_added_field(self):
        """ Gets bc/added data, handles var length barcodes plus added field"""

        # Should return shortest barcode, as this has the match to the added
        # field
        ids_bcs_added_field = {("AACTCGT", "1"): "s1",
                               ("AGCAGCACTTGT", "2"): "s2", ("AGCAG", "1"): "s3"}
        fasta_label = "123ABC region=1 length=255"
        fasta_seq = "AGCAGCACTTGTACAGATTAGACCGAGATTACGAGATTACCAAGAT"
        bc_lens = [12, 7, 5]
        all_bcs = ["AACTCGT", "AGCAGCACTTGT", "AGCAG"]
        barcode_type = "variable_length"
        max_bc_errors = 1.5
        disable_bc_correction = True
        added_demultiplex_field = 'region'
        curr_bc, corrected_bc, num_errors, added_field =\
            get_demultiplex_data(ids_bcs_added_field, fasta_label, fasta_seq,
                                 bc_lens, all_bcs, barcode_type, max_bc_errors, disable_bc_correction,
                                 added_demultiplex_field)

        expected_curr_bc = "AGCAG"
        expected_corrected_bc = "AGCAG"
        expected_num_errors = 0
        expected_added_field = '1'

        self.assertEqual(curr_bc, expected_curr_bc)
        self.assertEqual(corrected_bc, expected_corrected_bc)
        self.assertEqual(num_errors, expected_num_errors)
        self.assertEqual(added_field, expected_added_field)

    def test_write_qual_line_short_seq(self):
        """ Properly formats, writes qual score label/line for short seq """

        demultiplexed_qual_f = FakeOutFile()
        qual_seq = [25, 24, 22, 24, 24, 24, 25, 30, 23, 22, 22, 24, 25]
        label_line = "sample3_1 ABCD1234"
        keep_barcode = False
        bc_len = 4
        write_qual_line(demultiplexed_qual_f, qual_seq, label_line,
                        keep_barcode, bc_len)

        expected_data = '>sample3_1 ABCD1234\n24 24 25 30 23 22 22 24 25\n'

        self.assertEqual(demultiplexed_qual_f.data, expected_data)

    def test_write_qual_line_short_seq_retains_bc(self):
        """ Properly formats, writes qual score, retains BC segment """

        demultiplexed_qual_f = FakeOutFile()
        qual_seq = [25, 24, 22, 24, 24, 24, 25, 30, 23, 22, 22, 24, 25]
        label_line = "sample3_1 ABCD1234"
        keep_barcode = True
        bc_len = 4
        write_qual_line(demultiplexed_qual_f, qual_seq, label_line,
                        keep_barcode, bc_len)

        expected_data =\
            '>sample3_1 ABCD1234\n25 24 22 24 24 24 25 30 23 22 22 24 25\n'

        self.assertEqual(demultiplexed_qual_f.data, expected_data)

    def test_write_qual_line_long_seq(self):
        """ Properly formats, writes qual score label/line for long seq """

        demultiplexed_qual_f = FakeOutFile()
        qual_seq = [25, 24, 22, 24, 24, 24, 25, 30, 23, 22, 22, 24, 25,
                    14, 25, 27, 29, 30, 14, 10, 1, 23, 24, 27, 28, 30, 22, 24, 21,
                    24, 22, 21, 15, 17, 17, 15, 22, 13, 11, 10, 22, 24, 27, 28, 30,
                    14, 25, 27, 29, 30, 14, 10, 1, 23, 24, 27, 28, 30, 22, 24, 21,
                    14, 25, 27, 29, 30, 14, 10, 1, 23, 24, 27, 28, 30, 22, 24, 21]

        label_line = "sample3_1 ABCD1234"
        keep_barcode = False
        bc_len = 4
        write_qual_line(demultiplexed_qual_f, qual_seq, label_line,
                        keep_barcode, bc_len)

        expected_data = '>sample3_1 ABCD1234\n24 24 25 30 23 22 22 24 25 14 25 27 29 30 14 10 1 23 24 27 28 30 22 24 21 24 22 21 15 17 17 15 22 13 11 10 22 24 27 28 30 14 25 27 29 30 14 10 1 23 24 27 28 30 22 24 21 14 25 27\n29 30 14 10 1 23 24 27 28 30 22 24 21\n'
        self.assertEqual(demultiplexed_qual_f.data, expected_data)

    def test_write_qual_line_long_seq_retain_bc(self):
        """ Properly formats, writes qual long qual score with bc """

        demultiplexed_qual_f = FakeOutFile()
        qual_seq = [25, 24, 22, 24, 24, 24, 25, 30, 23, 22, 22, 24, 25,
                    14, 25, 27, 29, 30, 14, 10, 1, 23, 24, 27, 28, 30, 22, 24, 21,
                    24, 22, 21, 15, 17, 17, 15, 22, 13, 11, 10, 22, 24, 27, 28, 30,
                    14, 25, 27, 29, 30, 14, 10, 1, 23, 24, 27, 28, 30, 22, 24, 21,
                    14, 25, 27, 29, 30, 14, 10, 1, 23, 24, 27, 28, 30, 22, 24, 21]

        label_line = "sample3_1 ABCD1234"
        keep_barcode = True
        bc_len = 4
        write_qual_line(demultiplexed_qual_f, qual_seq, label_line,
                        keep_barcode, bc_len)

        expected_data = '>sample3_1 ABCD1234\n25 24 22 24 24 24 25 30 23 22 22 24 25 14 25 27 29 30 14 10 1 23 24 27 28 30 22 24 21 24 22 21 15 17 17 15 22 13 11 10 22 24 27 28 30 14 25 27 29 30 14 10 1 23 24 27 28 30 22 24\n21 14 25 27 29 30 14 10 1 23 24 27 28 30 22 24 21\n'
        self.assertEqual(demultiplexed_qual_f.data, expected_data)

    def test_write_fasta_line(self):
        """ Properly formats, write fasta label/line """

        demultiplexed_seqs_f = FakeOutFile()
        fasta_seq = "ACTAGACCTACAGGATACCATAGGACCAGATTTACA"
        label_line = "Sample1_213 ABCD1234"
        keep_barcode = False
        bc_len = 4
        write_fasta_line(demultiplexed_seqs_f, fasta_seq, label_line,
                         keep_barcode, bc_len)

        expected_data = ">Sample1_213 ABCD1234\nGACCTACAGGATACCATAGGACCAGATTTACA\n"
        self.assertEqual(demultiplexed_seqs_f.data, expected_data)

    def test_write_fasta_line_retain_bc(self):
        """ Properly formats, write fasta label/line, retains bc seqs """

        demultiplexed_seqs_f = FakeOutFile()
        fasta_seq = "ACTAGACCTACAGGATACCATAGGACCAGATTTACA"
        label_line = "Sample1_213 ABCD1234"
        keep_barcode = True
        bc_len = 4
        write_fasta_line(demultiplexed_seqs_f, fasta_seq, label_line,
                         keep_barcode, bc_len)

        expected_data = ">Sample1_213 ABCD1234\nACTAGACCTACAGGATACCATAGGACCAGATTTACA\n"
        self.assertEqual(demultiplexed_seqs_f.data, expected_data)

    def test_get_label_line(self):
        """ Properly formats fasta label """

        sample_id = "Sample1"
        fasta_label = "ABCD1234 region=1 length=254"
        bc = "AAAA"
        corrected_bc = "AAAT"
        num_errors = 1
        actual_label = get_label_line(sample_id, fasta_label, bc, corrected_bc,
                                      num_errors)

        expected_label = "Sample1 ABCD1234 orig_bc=AAAA new_bc=AAAT bc_diffs=1"
        self.assertEqual(actual_label, expected_label)

    def test_initialize_log_data(self):
        """ Properly initializes dict of zero counts for sample IDs """

        ids_bcs_added_field = {('AAAA', ''): 's1', ('TTTT', ''): 's2'}
        actual_log_data = initialize_log_data(ids_bcs_added_field)

        expected_log_data = {'TTTT,s2': 0, 'AAAA,s1': 0}

        self.assertEqual(actual_log_data, expected_log_data)

        # Handles added demultiplex field data
        ids_bcs_added_field = {('AAAA', '1'): 's1', ('TTTT', '2'): 's2'}
        actual_log_data = initialize_log_data(ids_bcs_added_field)

        expected_log_data = {'TTTT,2,s2': 0, 'AAAA,1,s1': 0}

        self.assertEqual(actual_log_data, expected_log_data)

    def test_get_output_ids_correct_bc(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', ''): 's1', ('TTTT', ''): 's2'}
        corrected_bc = 'AAAA'
        num_errors = 0
        added_field = None

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 's1_1'
        expected_log_id = 'AAAA,s1'
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_correct_bc_with_correction(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', ''): 's1', ('TTTT', ''): 's2'}
        corrected_bc = 'AAAA'
        num_errors = 1
        added_field = None

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 's1_1'
        expected_log_id = 'AAAA,s1'
        expected_bc_corrected_flag = 'corrected'

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_correct_bc_exceeds_max_errors(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', ''): 's1', ('TTTT', ''): 's2'}
        corrected_bc = 'AAAA'
        num_errors = 2
        added_field = None

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 'Unassigned_1'
        expected_log_id = ''
        expected_bc_corrected_flag = 'not_corrected'

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_bc_no_matches(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', ''): 's1', ('TTTT', ''): 's2'}
        corrected_bc = 'CCCC'
        num_errors = 0
        added_field = None

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 'Unassigned_1'
        expected_log_id = ''
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_no_bc_with_added_field(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('', '1'): 's1', ('', '2'): 's2'}
        corrected_bc = None
        num_errors = 0
        added_field = '1'

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 's1_1'
        expected_log_id = '1,s1'
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_no_bc_with_added_field_no_match(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('', '1'): 's1', ('', '2'): 's2'}
        corrected_bc = None
        num_errors = 0
        added_field = '3'

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 'Unassigned_1'
        expected_log_id = ''
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_correct_bc_added_field(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', '1'): 's1', ('TTTT', '2'): 's2'}
        corrected_bc = 'AAAA'
        num_errors = 0
        added_field = '1'

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 's1_1'
        expected_log_id = 'AAAA,1,s1'
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_correct_bc_added_field_no_match(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', '1'): 's1', ('TTTT', '2'): 's2'}
        corrected_bc = 'AAAA'
        num_errors = 0
        added_field = '3'

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 'Unassigned_1'
        expected_log_id = ''
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_get_output_ids_incorrect_bc_added_field_match(self):
        """ Properly formats SampleIDs/enumeration for writing seqs """

        ids_bcs_added_field = {('AAAA', '1'): 's1', ('TTTT', '2'): 's2'}
        corrected_bc = 'CCCC'
        num_errors = 0
        added_field = '1'

        sample_id, log_id, bc_corrected_flag =\
            get_output_ids(ids_bcs_added_field, corrected_bc, num_errors,
                           added_field)

        expected_sample_id = 'Unassigned_1'
        expected_log_id = ''
        expected_bc_corrected_flag = None

        self.assertEqual(sample_id, expected_sample_id)
        self.assertEqual(log_id, expected_log_id)
        self.assertEqual(bc_corrected_flag, expected_bc_corrected_flag)

    def test_assign_seqs_fasta_only(self):
        """ Properly iterates through, demultiplexes valid seq files """

        # Initial test for single fasta file alone.
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_no_errors]
        file_data['qual_files'] = []
        #file_data['mapping_file'] = self.valid_mapping_data_golay_upper
        file_data['demultiplexed_seqs_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = ">s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nCCAGATTACGAGATTA\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n"

        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 1, 'AACTCGTCGATG,s1': 1,
                             'AGCAGCACTTGT,s2': 1}
        expected_bc_freqs = {'AACTCGTCGATG': 1, 'AGCAGCACTTGT': 1,
                             'ACCGCAGAGTCA': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_fasta_plus_qual(self):
        """ Properly iterates through, demultiplexes valid seq files """

        # Handles single fasta and single qual
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_no_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = ">s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nCCAGATTACGAGATTA\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n"
        expected_demultiplexed_qual_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\n12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'

        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 1, 'AACTCGTCGATG,s1': 1,
                             'AGCAGCACTTGT,s2': 1}
        expected_bc_freqs = {'AACTCGTCGATG': 1, 'AGCAGCACTTGT': 1,
                             'ACCGCAGAGTCA': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_two_fastas(self):
        """ Properly iterates through, demultiplexes valid seq files """

        # Handles two fasta files alone
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_no_errors,
                                    self.valid_fasta_file_no_errors]
        file_data['qual_files'] = []
        #file_data['mapping_file'] = self.valid_mapping_data_golay_upper
        file_data['demultiplexed_seqs_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nCCAGATTACGAGATTA\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n>s1_4 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n>s3_5 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nCCAGATTACGAGATTA\n>s2_6 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 2, 'AACTCGTCGATG,s1': 2,
                             'AGCAGCACTTGT,s2': 2}
        expected_bc_freqs = {'AACTCGTCGATG': 2, 'AGCAGCACTTGT': 2,
                             'ACCGCAGAGTCA': 2}
        expected_seq_counts = 6
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_two_fastas_quals(self):
        """ Properly iterates through, demultiplexes valid seq files """

        # Handles single fasta and single qual
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_no_errors,
                                    self.valid_fasta_file_no_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors,
                                   self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nCCAGATTACGAGATTA\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n>s1_4 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n>s3_5 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nCCAGATTACGAGATTA\n>s2_6 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\n12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n>s1_4 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s3_5 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\n12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s2_6 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 2, 'AACTCGTCGATG,s1': 2,
                             'AGCAGCACTTGT,s2': 2}
        expected_bc_freqs = {'AACTCGTCGATG': 2, 'AGCAGCACTTGT': 2,
                             'ACCGCAGAGTCA': 2}
        expected_seq_counts = 6
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_error_correction(self):
        """ Properly iterates through, demultiplexes with error correction """

        # Handles single fasta and single qual with error correction
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_with_bc_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_1 ABCD0001 orig_bc=TACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=1\nCAGGACGAGACGAGGTT\n>s3_2 EFGH0002 orig_bc=GCCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=1\nCCAGATTACGAGATTA\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s1_1 ABCD0001 orig_bc=TACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=1\n30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s3_2 EFGH0002 orig_bc=GCCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=1\n12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 1, 'AACTCGTCGATG,s1': 1,
                             'AGCAGCACTTGT,s2': 1}
        expected_bc_freqs = {'TACTCGTCGATG': 1, 'GCCGCAGAGTCA': 1,
                             'AGCAGCACTTGT': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [2, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_exceeds_error_correction(self):
        """ Properly iterates through, demultiplexes with error correction """

        # Handles single fasta and single qual, errors exceed max
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_with_bc_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 0.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 0, 'AACTCGTCGATG,s1': 0,
                             'AGCAGCACTTGT,s2': 1}
        expected_bc_freqs = {'TACTCGTCGATG': 1, 'GCCGCAGAGTCA': 1,
                             'AGCAGCACTTGT': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 2]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_exceeds_error_correction_unassigned(self):
        """ Properly iterates through, demultiplexes with error correction """

        # Handles single fasta and single qual, disabled bc correction,
        # writes unassigned sequences, retains barcodes, starts enumeration
        # at 1000, generic 12 base pair barcode type.
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_with_bc_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()
        file_data['unassigned_seqs_f'] = FakeOutFile()
        file_data['unassigned_qual_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', ''): 's1',
                               ('AGCAGCACTTGT', ''): 's2', ('ACCGCAGAGTCA', ''): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = True
        barcode_type = 12
        max_bc_errors = 1.5
        start_index = 1000
        write_unassigned_reads = True
        disable_bc_correction = True
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s2_1002 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nAGCAGCACTTGTGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s2_1002 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)
        expected_unassigned_fasta_seq = '>Unassigned_1000 ABCD0001 orig_bc=TACTCGTCGATG new_bc=None bc_diffs=0\nTACTCGTCGATGCAGGACGAGACGAGGTT\n>Unassigned_1001 EFGH0002 orig_bc=GCCGCAGAGTCA new_bc=None bc_diffs=0\nGCCGCAGAGTCACCAGATTACGAGATTA\n'
        expected_unassigned_qual_seq = '>Unassigned_1000 ABCD0001 orig_bc=TACTCGTCGATG new_bc=None bc_diffs=0\n29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>Unassigned_1001 EFGH0002 orig_bc=GCCGCAGAGTCA new_bc=None bc_diffs=0\n13 22 15 12 10 14 23 13 25 22 15 20 12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n'
        self.assertEqual(file_data['unassigned_seqs_f'].data,
                         expected_unassigned_fasta_seq)
        self.assertEqual(file_data['unassigned_qual_f'].data,
                         expected_unassigned_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 0, 'AACTCGTCGATG,s1': 0,
                             'AGCAGCACTTGT,s2': 1}
        expected_bc_freqs = {'TACTCGTCGATG': 1, 'GCCGCAGAGTCA': 1,
                             'AGCAGCACTTGT': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_assign_seqs_fasta_qual_added_demultiplex_field(self):
        """ Properly iterates through, demultiplexes valid seq files """

        # Handles single fasta and single qual
        file_data = {}
        file_data['fasta_files'] = [self.valid_fasta_file_no_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()

        ids_bcs_added_field = {('AACTCGTCGATG', '1'): 's1',
                               ('AGCAGCACTTGT', '2'): 's2', ('ACCGCAGAGTCA', '3'): 's3'}
        bc_lens = [12]
        all_bcs = ['AACTCGTCGATG', 'AGCAGCACTTGT', 'ACCGCAGAGTCA']
        keep_barcode = False
        barcode_type = "golay_12"
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = 'Added_Demultiplex'

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            assign_seqs(file_data, ids_bcs_added_field, bc_lens, all_bcs,
                        keep_barcode, barcode_type, max_bc_errors, start_index,
                        write_unassigned_reads, disable_bc_correction,
                        added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nCAGGACGAGACGAGGTT\n'
        expected_demultiplexed_qual_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n'

        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,3,s3': 0, 'AGCAGCACTTGT,2,s2': 0,
                             'AACTCGTCGATG,1,s1': 1}
        expected_bc_freqs = {'AACTCGTCGATG': 1, 'AGCAGCACTTGT': 1,
                             'ACCGCAGAGTCA': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_demultiplex_sequences_default_settings(self):
        """ Overall functionality test """

        file_data = {}
        file_data['mapping_file'] = self.valid_mapping_data_golay_upper
        file_data['fasta_files'] = [self.valid_fasta_file_no_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()

        keep_barcode = False,
        barcode_type = 'golay_12'
        max_bc_errors = 1.5
        start_index = 1
        write_unassigned_reads = False
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            demultiplex_sequences(file_data, keep_barcode, barcode_type,
                                  max_bc_errors, start_index, write_unassigned_reads,
                                  disable_bc_correction, added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\nAACTCGTCGATGCAGGACGAGACGAGGTT\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\nACCGCAGAGTCACCAGATTACGAGATTA\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\nAGCAGCACTTGTGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s1_1 ABCD0001 orig_bc=AACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=0\n29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s3_2 EFGH0002 orig_bc=ACCGCAGAGTCA new_bc=ACCGCAGAGTCA bc_diffs=0\n13 22 15 12 10 14 23 13 25 22 15 20 12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s2_3 IJKL0003 orig_bc=AGCAGCACTTGT new_bc=AGCAGCACTTGT bc_diffs=0\n30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'ACCGCAGAGTCA,s3': 1, 'AACTCGTCGATG,s1': 1,
                             'AGCAGCACTTGT,s2': 1}
        expected_bc_freqs = {'AACTCGTCGATG': 1, 'AGCAGCACTTGT': 1,
                             'ACCGCAGAGTCA': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_demultiplex_sequences_alternate_settings(self):
        """ Overall functionality test with alternate settings"""

        file_data = {}
        file_data['mapping_file'] = self.valid_mapping_data_golay_upper
        file_data['fasta_files'] = [self.valid_fasta_file_with_bc_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()
        file_data['unassigned_seqs_f'] = FakeOutFile()
        file_data['unassigned_qual_f'] = FakeOutFile()

        keep_barcode = True,
        barcode_type = 12
        max_bc_errors = 1
        start_index = 500
        write_unassigned_reads = True
        disable_bc_correction = False
        added_demultiplex_field = 'Added_Demultiplex'

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            demultiplex_sequences(file_data, keep_barcode, barcode_type,
                                  max_bc_errors, start_index, write_unassigned_reads,
                                  disable_bc_correction, added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_500 ABCD0001 orig_bc=TACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=1\nTACTCGTCGATGCAGGACGAGACGAGGTT\n'
        expected_demultiplexed_qual_seq = '>s1_500 ABCD0001 orig_bc=TACTCGTCGATG new_bc=AACTCGTCGATG bc_diffs=1\n29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'AGCAGCACTTGT,1,s2': 0, 'ACCGCAGAGTCA,1,s3': 0,
                             'AACTCGTCGATG,1,s1': 1}
        expected_bc_freqs = {'TACTCGTCGATG': 1, 'GCCGCAGAGTCA': 1,
                             'AGCAGCACTTGT': 1}
        expected_seq_counts = 3
        expected_corrected_bc_count = [1, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_demultiplex_sequences_single_sample_no_demultiplex_data(self):
        """ Overall functionality with a single sample, no demultiplex data"""

        file_data = {}
        file_data['mapping_file'] = \
            self.valid_mapping_data_no_bcs_no_added_demultiplex
        file_data['fasta_files'] = [self.valid_fasta_file_with_bc_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()
        file_data['unassigned_seqs_f'] = FakeOutFile()
        file_data['unassigned_qual_f'] = FakeOutFile()

        keep_barcode = False,
        barcode_type = 0
        max_bc_errors = 1
        start_index = 500
        write_unassigned_reads = True
        disable_bc_correction = False
        added_demultiplex_field = None

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            demultiplex_sequences(file_data, keep_barcode, barcode_type,
                                  max_bc_errors, start_index, write_unassigned_reads,
                                  disable_bc_correction, added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_500 ABCD0001 orig_bc= new_bc= bc_diffs=0\nTACTCGTCGATGCAGGACGAGACGAGGTT\n>s1_501 EFGH0002 orig_bc= new_bc= bc_diffs=0\nGCCGCAGAGTCACCAGATTACGAGATTA\n>s1_502 IJKL0003 orig_bc= new_bc= bc_diffs=0\nAGCAGCACTTGTGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s1_500 ABCD0001 orig_bc= new_bc= bc_diffs=0\n29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s1_501 EFGH0002 orig_bc= new_bc= bc_diffs=0\n13 22 15 12 10 14 23 13 25 22 15 20 12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s1_502 IJKL0003 orig_bc= new_bc= bc_diffs=0\n30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'s1': 3}
        expected_bc_freqs = {'': 3}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

    def test_demultiplex_sequences_added_demultiplex_only(self):
        """ Overall functionality test with alternate settings"""

        file_data = {}
        file_data['mapping_file'] =\
            self.valid_mapping_data_no_bcs_added_demultiplex
        file_data['fasta_files'] = [self.valid_fasta_file_with_bc_errors]
        file_data['qual_files'] = [self.valid_qual_file_no_errors]
        file_data['demultiplexed_seqs_f'] = FakeOutFile()
        file_data['demultiplexed_qual_f'] = FakeOutFile()
        file_data['unassigned_seqs_f'] = FakeOutFile()
        file_data['unassigned_qual_f'] = FakeOutFile()

        keep_barcode = True,
        barcode_type = 0
        max_bc_errors = 1
        start_index = 500
        write_unassigned_reads = True
        disable_bc_correction = False
        added_demultiplex_field = 'Added_Demultiplex'

        log_data, bc_freqs, seq_counts, corrected_bc_count =\
            demultiplex_sequences(file_data, keep_barcode, barcode_type,
                                  max_bc_errors, start_index, write_unassigned_reads,
                                  disable_bc_correction, added_demultiplex_field)

        expected_demultiplexed_fasta_seq = '>s1_500 ABCD0001 orig_bc= new_bc= bc_diffs=0\nTACTCGTCGATGCAGGACGAGACGAGGTT\n>s2_501 EFGH0002 orig_bc= new_bc= bc_diffs=0\nGCCGCAGAGTCACCAGATTACGAGATTA\n>s3_502 IJKL0003 orig_bc= new_bc= bc_diffs=0\nAGCAGCACTTGTGACCGATTACGATAACG\n'
        expected_demultiplexed_qual_seq = '>s1_500 ABCD0001 orig_bc= new_bc= bc_diffs=0\n29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24\n>s2_501 EFGH0002 orig_bc= new_bc= bc_diffs=0\n13 22 15 12 10 14 23 13 25 22 15 20 12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13\n>s3_502 IJKL0003 orig_bc= new_bc= bc_diffs=0\n30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17\n'
        self.assertEqual(file_data['demultiplexed_seqs_f'].data,
                         expected_demultiplexed_fasta_seq)
        self.assertEqual(file_data['demultiplexed_qual_f'].data,
                         expected_demultiplexed_qual_seq)

        expected_log_data = {'2,s2': 1, '3,s3': 1, '1,s1': 1}
        expected_bc_freqs = {'': 3}
        expected_seq_counts = 3
        expected_corrected_bc_count = [0, 0]

        self.assertEqual(log_data, expected_log_data)
        self.assertEqual(bc_freqs, expected_bc_freqs)
        self.assertEqual(seq_counts, expected_seq_counts)
        self.assertEqual(corrected_bc_count, expected_corrected_bc_count)

# Large test data strings at the end for better readability


valid_mapping_data_golay =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\ts1_description',
     's2\tagcAGCACTTGT\tATTCGATART\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']

valid_mapping_data_golay_upper =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tAdded_Demultiplex\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\t1\ts1_description',
     's2\tAGCAGCACTTGT\tATTCGATART\t1\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\t1\ts3_description']

# BarcodeSequence header not listed correctly
invalid_mapping_data_golay_header =\
    ['#SampleID\tBarcode\tLinkerPrimerSequence\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\ts1_description',
     's2\tAGCAGCACTTGT\tATTCGATART\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']
# non-DNA characters in barcode
invalid_mapping_data_golay_dna =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\tAACTCGTCGAXG\tATTCGATART\ts1_description',
     's2\tAGCAGCACTTGT\tATTCGATART\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']
# duplicate barcodes
invalid_mapping_data_golay_dup_bcs =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\tACCGCAGAGTCA\tATTCGATART\ts1_description',
     's2\tAGCAGCACTTGT\tATTCGATART\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']
# duplicate SampleIDs
invalid_mapping_data_golay_dup_sids =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\ts1_description',
     's2\tAGCAGCACTTGT\tATTCGATART\ts2_description',
     's2\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']
# Added demultiplex has duplicates
invalid_mapping_data_no_bcs_added_demultiplex =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tAdded_Demultiplex\tDescription',
     's1\t\tATTCGATART\t1\ts1_description',
     's2\t\tATTCGATART\t2\ts2_description',
     's3\t\tYATGCTGCCTCCCGTAGGAGT\t1\ts3_description']
# barcodes plus added demultiplex still has duplicates
invalid_mapping_data_bcs_added_demultiplex =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tAdded_Demultiplex\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\t1\ts1_description',
     's2\tAACTCCCCGATG\tATTCGATART\t2\ts2_description',
     's3\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\t1\ts3_description']
# Missing barcode
invalid_mapping_data_golay_missing_bc =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\ts1_description',
     's2\t\tATTCGATART\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']

valid_mapping_data_var_len_bcs =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\tAACTCGTCGATG\tATTCGATART\ts1_description',
     's2\tCACTTGT\tATTCGATART\ts2_description',
     's3\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\ts3_description']

valid_mapping_data_no_bcs_added_demultiplex =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tAdded_Demultiplex\tDescription',
     's1\t\tATTCGATART\t1\ts1_description',
     's2\t\tATTCGATART\t2\ts2_description',
     's3\t\tYATGCTGCCTCCCGTAGGAGT\t3\ts3_description']

valid_mapping_data_bcs_and_added_demultiplex =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tAdded_Demultiplex\tDescription',
     's1\tAAAA\tATTCGATART\t1\ts1_description',
     's2\tAAAA\tATTCGATART\t2\ts2_description',
     's3\tCCCC\tYATGCTGCCTCCCGTAGGAGT\t1\ts3_description']

# Single sample with no demultiplexing data
valid_mapping_data_no_bcs_no_added_demultiplex =\
    ['#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription',
     's1\t\t\ts1_description']

expected_formatted_log_data_no_qual = ['demultiplex_fasta.py log data\n',
                                       'Metadata mapping file:\ttest_mapping.txt',
                                       'Input FASTA file(s):\tfasta1.fna,fasta2.fna',
                                       'Total sequences in input files:\t32',
                                       'Retain barcode:\tFalse', 'Barcode type:\tgolay_12',
                                       'Max barcode error/mismatches allowed:\t1.5',
                                       'Starting sequence identifier:\t1',
                                       'Write unassigned reads:\tFalse',
                                       'Disable barcode correction:\tFalse',
                                       'Added demultiplex field:\tNone',
                                       'Save barcode frequencies:\tFalse\n',
                                       'Barcodes corrected/not corrected:\t7/15',
                                       'Number of samples in mapping file:\t3',
                                       'Sample count min/max/mean:\t0 / 20 / 10.00',
                                       'Sample\tSequence Count\tBarcode/Added Demultiplex',
                                       'TTTT\t20\t', 'AAAA\t10\t', 'CCCC\t0\t', 'Seqs written\t30',
                                       'Percent of input seqs written\t0.94']

expected_formatted_log_data_with_qual = ['demultiplex_fasta.py log data\n',
                                         'Metadata mapping file:\ttest_mapping.txt',
                                         'Input FASTA file(s):\tfasta1.fna',
                                         'Input QUAL file(s):\tqual1.qual',
                                         'Total sequences in input files:\t32',
                                         'Retain barcode:\tTrue', 'Barcode type:\t8',
                                         'Max barcode error/mismatches allowed:\t0',
                                         'Starting sequence identifier:\t1000', 'Write unassigned reads:\tTrue',
                                         'Disable barcode correction:\tTrue', 'Added demultiplex field:\trun_prefix',
                                         'Save barcode frequencies:\tTrue\n', 'Barcodes corrected/not corrected:\t7/15',
                                         'Number of samples in mapping file:\t3',
                                         'Sample count min/max/mean:\t0 / 20 / 10.00',
                                         'Sample\tSequence Count\tBarcode/Added Demultiplex', 'TTTT\t20\t',
                                         'AAAA\t10\t', 'CCCC\t0\t', 'Seqs written\t30',
                                         'Percent of input seqs written\t0.94']


valid_fasta_file_no_errors = ['>ABCD0001 Added_Demultiplex=1 length=254',
                              'AACTCGTCGATGCAGGACGAGACGAGGTT',
                              '>EFGH0002 Added_Demultiplex=2 length=254',
                              'ACCGCAGAGTCACCAGATTACGAGATTA',
                              '>IJKL0003 Added_Demultiplex=3 length=255',
                              'AGCAGCACTTGTGACCGATTACGATAACG']

valid_qual_file_no_errors = ['>ABCD0001 Added_Demultiplex=1 length=254',
                             '29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24',
                             '>EFGH0002 Added_Demultiplex=2 length=254',
                             '13 22 15 12 10 14 23 13 25 22 15 20 12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13',
                             '>IJKL0003 Added_Demultiplex=3 length=255',
                             '30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17']

valid_fasta_file_with_bc_errors = ['>ABCD0001 Added_Demultiplex=1 length=254',
                                   'TACTCGTCGATGCAGGACGAGACGAGGTT',
                                   '>EFGH0002 Added_Demultiplex=2 length=254',
                                   'GCCGCAGAGTCACCAGATTACGAGATTA',
                                   '>IJKL0003 Added_Demultiplex=3 length=255',
                                   'AGCAGCACTTGTGACCGATTACGATAACG']

sample_correct_mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	ReversePrimer	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356"""

sample_fasta_file = """>ABCD0001 Added_Demultiplex=1 length=254
AGCACGAGCCTACAGGACGAGACGAGGTT
>EFGH0002 Added_Demultiplex=2 length=254
AACTCGTCGATGCCAGATTACGAGATTA
>IJKL0003 Added_Demultiplex=3 length=255
ACAGACCACTCAGACCGATTACGATAACG
>nomatchingBC Added_Demultiplex=3 length=255
ACCCCCCACTCAGACCGATTACGATAACG"""

sample_qual_file = """>ABCD0001 Added_Demultiplex=1 length=254
29 13 24 14 10 14 16 13 30 10 13 11 30 26 11 11 29 20 19 16 24 17 29 28 11 27 14 24 24
>EFGH0002 Added_Demultiplex=2 length=254
13 22 15 12 10 14 23 13 25 22 15 20 12 14 27 23 22 19 24 18 19 20 28 10 17 14 17 13
>IJKL0003 Added_Demultiplex=3 length=255
30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17
>nomatchingBC Added_Demultiplex=3 length=255
30 27 11 16 30 19 13 19 16 15 24 12 10 20 16 20 25 27 22 28 16 22 16 18 12 13 16 25 17"""

if __name__ == '__main__':
    main()
