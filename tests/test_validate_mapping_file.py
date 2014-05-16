#!/usr/bin/env python

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["Rob Knight", "William Walters"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"

from os import close
from os.path import isdir, isfile, exists, join, basename
from tempfile import mkdtemp
from shutil import rmtree
from collections import defaultdict
from tempfile import mkstemp, mkdtemp

from unittest import TestCase, main
from skbio.util.misc import remove_files
from skbio.util.misc import create_dir

from qiime.check_id_map import (check_mapping_file, process_id_map,
                                check_data_fields, check_fields_past_bounds, check_chars_data_fields,
                                check_dna_chars_primers, check_dna_chars_bcs, check_bcs_lengths,
                                check_bc_duplicates, check_fixed_len_bcs_dups, check_variable_len_bcs_dups,
                                check_added_demultiplex_dups, check_sampleid_duplicates, check_header,
                                check_header_dups, check_header_chars, check_header_required_fields,
                                correct_mapping_data)


class CheckIdMapTests(TestCase):

    def setUp(self):
        # create the temporary input files that will be used

        self._files_to_remove = []

        # Input data
        self.sample_correct_mapping_data = sample_correct_mapping_data
        self.sample_errors_mapping_data = sample_errors_mapping_data
        self.sample_warnings_mapping_data = sample_warnings_mapping_data
        self.sample_errors_warnings_mapping_data =\
            sample_errors_warnings_mapping_data
        self.empty_fields_mapping_data = empty_fields_mapping_data

        fd, self.correct_mapping_fp = mkstemp(
            prefix='correct_mapping_',
            suffix='.txt')
        close(fd)
        map_file = open(self.correct_mapping_fp, 'w')
        map_file.write(self.sample_correct_mapping_data)
        map_file.close()

        fd, self.errors_mapping_fp = mkstemp(
            prefix='errors_mapping_',
            suffix='.txt')
        close(fd)
        map_file = open(self.errors_mapping_fp, 'w')
        map_file.write(self.sample_errors_mapping_data)
        map_file.close()

        fd, self.empty_fields_fp = mkstemp(
            prefix='empty_fields_',
            suffix='.txt')
        close(fd)
        map_file = open(self.empty_fields_fp, 'w')
        map_file.write(self.empty_fields_mapping_data)
        map_file.close()

        fd, self.warnings_mapping_fp = mkstemp(
            prefix='warnings_mapping_',
            suffix='.txt')
        close(fd)
        map_file = open(self.warnings_mapping_fp, 'w')
        map_file.write(self.sample_warnings_mapping_data)
        map_file.close()

        fd, self.errors_warnings_mapping_fp = mkstemp(
            prefix='errors_warnings_mapping_',
            suffix='.txt')
        close(fd)
        map_file = open(self.errors_warnings_mapping_fp, 'w')
        map_file.write(self.sample_errors_warnings_mapping_data)
        map_file.close()

        # Output data
        self.expected_html_data_correct_input = expected_html_data_correct_input
        self.expected_corrected_data_correct_input =\
            expected_corrected_data_correct_input
        self.expected_log_data_correct_input = expected_log_data_correct_input

        self.expected_html_errors_output = expected_html_errors_output
        self.expected_data_errors_corrected_output =\
            expected_data_errors_corrected_output
        self.expected_data_log_errors_output = expected_data_log_errors_output
        self.expected_html_errors_suppressed_bcs =\
            expected_html_errors_suppressed_bcs
        self.expected_output_log_errors_bcs_suppressed =\
            expected_output_log_errors_bcs_suppressed

        self.expected_html_output_warnings = expected_html_output_warnings
        self.expected_corrected_warnings_output =\
            expected_corrected_warnings_output
        self.expected_log_warnings_output =\
            expected_log_warnings_output

        self.expected_html_errors_warnings_output =\
            expected_html_errors_warnings_output
        self.expected_corrected_data_errors_warnings =\
            expected_corrected_data_errors_warnings
        self.expected_log_errors_warnings_output =\
            expected_log_errors_warnings_output

        self.output_dir = mkdtemp()
        self.output_dir += '/'

        create_dir(self.output_dir)

        self._files_to_remove =\
            [self.correct_mapping_fp, self.errors_mapping_fp,
             self.warnings_mapping_fp, self.errors_warnings_mapping_fp,
             self.empty_fields_fp]

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)

    def test_check_mapping_file_correct_file(self):
        """ Gives proper files for valid mapping files """

        # Use valid data, default parameters
        check_mapping_file(mapping_fp=self.correct_mapping_fp,
                           output_dir=self.output_dir,
                           verbose=False)

        # Check existence of expected output files
        output_html_fp = join(self.output_dir,
                              basename(self.correct_mapping_fp).replace('.txt', '.html'))
        output_corrected_fp =\
            join(self.output_dir,
                 basename(self.correct_mapping_fp).replace('.txt', '_corrected.txt'))
        output_log_fp =\
            join(self.output_dir,
                 basename(self.correct_mapping_fp).replace('.txt', '.log'))
        overlib_js_fp = join(self.output_dir, 'overlib.js')

        self.assertTrue(exists(output_html_fp))
        self.assertTrue(exists(output_corrected_fp))
        self.assertTrue(exists(output_log_fp))
        self.assertTrue(exists(overlib_js_fp))

        # Check output data for expected results

        html_data = "".join([line for line in open(output_html_fp, "U")])
        corrected_data =\
            "".join([line for line in open(output_corrected_fp, "U")])
        log_data = "".join([line for line in open(output_log_fp, "U")])

        self.assertEqual(html_data, self.expected_html_data_correct_input)
        self.assertEqual(corrected_data,
                         self.expected_corrected_data_correct_input)
        self.assertEqual(log_data, self.expected_log_data_correct_input)

        # With additional parameters added should not change results using
        # same valid input data
        check_mapping_file(mapping_fp=self.correct_mapping_fp,
                           output_dir=self.output_dir,
                           has_barcodes=True,
                           char_replace="A",
                           verbose=False,
                           variable_len_barcodes=True,
                           disable_primer_check=True,
                           added_demultiplex_field=None)

        # Check existence of expected output files
        output_html_fp = join(self.output_dir,
                              basename(self.correct_mapping_fp).replace('.txt', '.html'))
        output_corrected_fp =\
            join(self.output_dir,
                 basename(self.correct_mapping_fp).replace('.txt', '_corrected.txt'))
        output_log_fp =\
            join(self.output_dir,
                 basename(self.correct_mapping_fp).replace('.txt', '.log'))
        overlib_js_fp = join(self.output_dir, 'overlib.js')

        self.assertTrue(exists(output_html_fp))
        self.assertTrue(exists(output_corrected_fp))
        self.assertTrue(exists(output_log_fp))
        self.assertTrue(exists(overlib_js_fp))

        # Check output data for expected results

        html_data = "".join([line for line in open(output_html_fp, "U")])
        corrected_data =\
            "".join([line for line in open(output_corrected_fp, "U")])
        log_data = "".join([line for line in open(output_log_fp, "U")])

        self.assertEqual(html_data, self.expected_html_data_correct_input)
        self.assertEqual(corrected_data,
                         self.expected_corrected_data_correct_input)
        self.assertEqual(log_data, self.expected_log_data_correct_input)

    def test_check_mapping_file_errors(self):
        """ Gives proper files for errors in mapping files """

        # Use data with errors, default parameters
        check_mapping_file(mapping_fp=self.errors_mapping_fp,
                           output_dir=self.output_dir,
                           verbose=False)

        # Check existence of expected output files
        output_html_fp = join(self.output_dir,
                              basename(self.errors_mapping_fp).replace('.txt', '.html'))
        output_corrected_fp =\
            join(self.output_dir,
                 basename(self.errors_mapping_fp).replace('.txt', '_corrected.txt'))
        output_log_fp =\
            join(self.output_dir,
                 basename(self.errors_mapping_fp).replace('.txt', '.log'))
        overlib_js_fp = join(self.output_dir, 'overlib.js')

        self.assertTrue(exists(output_html_fp))
        self.assertTrue(exists(output_corrected_fp))
        self.assertTrue(exists(output_log_fp))
        self.assertTrue(exists(overlib_js_fp))

        # Check output data for expected results

        html_data = "".join([line for line in open(output_html_fp, "U")])
        corrected_data =\
            "".join([line for line in open(output_corrected_fp, "U")])
        log_data = "".join([line for line in open(output_log_fp, "U")])

        self.assertEqual(html_data, self.expected_html_errors_output)
        self.assertEqual(corrected_data,
                         self.expected_data_errors_corrected_output)
        self.assertEqual(log_data, self.expected_data_log_errors_output)

    def test_check_mapping_file_errors_suppressed_bcs(self):
        """ Should suppress errors about barcodes if bcs turned off """

        # Should not flag bcs for errors with invalid characters
        check_mapping_file(mapping_fp=self.errors_mapping_fp,
                           output_dir=self.output_dir,
                           has_barcodes=False,
                           char_replace="A",
                           verbose=False,
                           variable_len_barcodes=True,
                           disable_primer_check=True,
                           added_demultiplex_field=None)

        # Check existence of expected output files
        output_html_fp = join(self.output_dir,
                              basename(self.errors_mapping_fp).replace('.txt', '.html'))
        output_corrected_fp =\
            join(self.output_dir,
                 basename(self.errors_mapping_fp).replace('.txt', '_corrected.txt'))
        output_log_fp =\
            join(self.output_dir,
                 basename(self.errors_mapping_fp).replace('.txt', '.log'))
        overlib_js_fp = join(self.output_dir, 'overlib.js')

        self.assertTrue(exists(output_html_fp))
        self.assertTrue(exists(output_corrected_fp))
        self.assertTrue(exists(output_log_fp))
        self.assertTrue(exists(overlib_js_fp))

        # Check output data for expected results

        html_data = "".join([line for line in open(output_html_fp, "U")])
        corrected_data =\
            "".join([line for line in open(output_corrected_fp, "U")])
        log_data = "".join([line for line in open(output_log_fp, "U")])

        self.assertEqual(html_data, self.expected_html_errors_suppressed_bcs)
        self.assertEqual(corrected_data,
                         self.expected_data_errors_corrected_output)
        self.assertEqual(log_data,
                         self.expected_output_log_errors_bcs_suppressed)

    def test_check_mapping_file_warnings(self):
        """ Gives proper files for warnings in mapping files """

        check_mapping_file(mapping_fp=self.warnings_mapping_fp,
                           output_dir=self.output_dir,
                           verbose=False)

        # Check existence of expected output files
        output_html_fp = join(self.output_dir,
                              basename(self.warnings_mapping_fp).replace('.txt', '.html'))
        output_corrected_fp =\
            join(self.output_dir,
                 basename(self.warnings_mapping_fp).replace('.txt', '_corrected.txt'))
        output_log_fp =\
            join(self.output_dir,
                 basename(self.warnings_mapping_fp).replace('.txt', '.log'))
        overlib_js_fp = join(self.output_dir, 'overlib.js')

        self.assertTrue(exists(output_html_fp))
        self.assertTrue(exists(output_corrected_fp))
        self.assertTrue(exists(output_log_fp))
        self.assertTrue(exists(overlib_js_fp))

        # Check output data for expected results

        html_data = "".join([line for line in open(output_html_fp, "U")])
        corrected_data =\
            "".join([line for line in open(output_corrected_fp, "U")])
        log_data = "".join([line for line in open(output_log_fp, "U")])

        self.assertEqual(html_data, self.expected_html_output_warnings)
        self.assertEqual(corrected_data,
                         self.expected_corrected_warnings_output)
        self.assertEqual(log_data, self.expected_log_warnings_output)

    def test_check_mapping_file_multiple_problems(self):
        """ Gives proper files for combinations of errors/warnings """

        check_mapping_file(mapping_fp=self.errors_warnings_mapping_fp,
                           output_dir=self.output_dir,
                           added_demultiplex_field="DoesNotExist",
                           verbose=False)

        # Check existence of expected output files
        output_html_fp = join(self.output_dir,
                              basename(self.errors_warnings_mapping_fp).replace('.txt', '.html'))
        output_corrected_fp =\
            join(self.output_dir,
                 basename(self.errors_warnings_mapping_fp).replace('.txt',
                                                                   '_corrected.txt'))
        output_log_fp =\
            join(self.output_dir,
                 basename(self.errors_warnings_mapping_fp).replace('.txt', '.log'))
        overlib_js_fp = join(self.output_dir, 'overlib.js')

        self.assertTrue(exists(output_html_fp))
        self.assertTrue(exists(output_corrected_fp))
        self.assertTrue(exists(output_log_fp))
        self.assertTrue(exists(overlib_js_fp))

        # Check output data for expected results

        html_data = "".join([line for line in open(output_html_fp, "U")])
        corrected_data =\
            "".join([line for line in open(output_corrected_fp, "U")])
        log_data = "".join([line for line in open(output_log_fp, "U")])

        self.assertEqual(html_data, self.expected_html_errors_warnings_output)
        self.assertEqual(corrected_data,
                         self.expected_corrected_data_errors_warnings)
        self.assertEqual(log_data, self.expected_log_errors_warnings_output)

    def test_process_id_map_correct_data(self):
        """ Returns expected results for correct mapping data """

        header, mapping_data, comments, errors, warnings =\
            process_id_map(self.correct_mapping_fp)

        expected_header = [
            'SampleID',
            'BarcodeSequence',
            'LinkerPrimerSequence',
            'Treatment',
            'ReversePrimer',
            'Description']
        expected_mapping_data = [['PC.354',
                                  'AGCACGAGCCTA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._354'],
                                 ['PC.355',
                                  'AACTCGTCGATG',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._355'],
                                 ['PC.356',
                                  'ACAGACCACTCA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._356']]
        expected_comments = [
            'Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).']
        expected_errors = []
        expected_warnings = []

        self.assertEqual(header, expected_header)
        self.assertEqual(mapping_data, expected_mapping_data)
        self.assertEqual(comments, expected_comments)
        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_process_id_map_errors(self):
        """ Returns expected results for mapping data with errors """

        header, mapping_data, comments, errors, warnings =\
            process_id_map(self.errors_mapping_fp)

        expected_header = [
            'SampleID',
            'BarcodeSequence',
            'LinkerPrimerSequence',
            'Treatment',
            'ReversePrimer',
            'NotDescription']
        expected_mapping_data = [['PC.355',
                                  'AGCACGAGCCxTA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._354'],
                                 ['PC.355',
                                  'AACTCGTCGATG',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._355'],
                                 ['PC.356',
                                  'ACAGACCACTCA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._356']]
        expected_comments = [
            'Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).']
        expected_errors = [
            'Found header field NotDescription, last field should be Description\t0,5',
            'Invalid DNA sequence detected: AGCACGAGCCxTA\t1,1',
            'Duplicate SampleID PC.355 found.\t1,0',
            'Duplicate SampleID PC.355 found.\t2,0']
        expected_warnings = [
            'Barcode AGCACGAGCCxTA differs than length 12\t1,1']

        self.assertEqual(header, expected_header)
        self.assertEqual(mapping_data, expected_mapping_data)
        self.assertEqual(comments, expected_comments)
        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_process_id_map_empty_data_fields(self):
        """ Returns expected results for mapping data with missing fields """

        header, mapping_data, comments, errors, warnings =\
            process_id_map(self.empty_fields_fp)

        expected_header = [
            'SampleID',
            'BarcodeSequence',
            'LinkerPrimerSequence',
            'Treatment',
            'ReversePrimer',
            'Description']
        expected_mapping_data = [['PC.354',
                                  'AGCACGAGCCTA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._354'],
                                 ['PC.355',
                                  'AACTCGTCGATG',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  '',
                                  ''],
                                 ['PC.356',
                                  'ACAGACCACTCA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._356']]
        expected_comments = [
            'Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).']
        expected_errors = ['Missing expected DNA sequence\t2,4']
        expected_warnings = [
            'Empty data field  found\t2,4',
            'Empty data field  found\t2,5']

        self.assertEqual(header, expected_header)
        self.assertEqual(mapping_data, expected_mapping_data)
        self.assertEqual(comments, expected_comments)
        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_process_id_map_warnings(self):
        """ Returns expected results for mapping data with warnings """

        header, mapping_data, comments, errors, warnings =\
            process_id_map(self.warnings_mapping_fp)

        expected_header = [
            'SampleID',
            'BarcodeSequence',
            'LinkerPrimerSequence',
            'Treatm-ent',
            'ReversePrimer',
            'Description']
        expected_mapping_data = [['PC.354',
                                  'AGCACGAGCCTA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._354'],
                                 ['PC_355',
                                  'AACTCGTCGATG',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Co&ntrol',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._355',
                                  'OutOfBounds'],
                                 ['PC.356',
                                  'ACAGACCACTCA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._356']]
        expected_comments = [
            'Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).']
        expected_errors = []
        expected_warnings = [
            'Found invalid character in Treatm-ent header field.\t0,3',
            'Invalid characters found in PC_355\t2,0',
            'Invalid characters found in Co&ntrol\t2,3',
            'Data field OutOfBounds found after Description column\t2,6']

        self.assertEqual(header, expected_header)
        self.assertEqual(mapping_data, expected_mapping_data)
        self.assertEqual(comments, expected_comments)
        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_process_id_map_multiple_problems(self):
        """ Returns expected results for combinations of problems """

        header, mapping_data, comments, errors, warnings =\
            process_id_map(self.errors_warnings_mapping_fp)

        expected_header = [
            'SampleID',
            'BarcodeSequence',
            'LinkerPrimerSequence',
            'Treatment',
            'Treatment',
            'Description']
        expected_mapping_data = [['PC.354',
                                  'AGCACGAGCCTA',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Cont^^rol',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._354'],
                                 ['PC-355',
                                  'AACTCGTCGATGN',
                                  'YATGCTGCCTCCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._355',
                                  'outofbounds'],
                                 ['PC.356',
                                  'ACAGACCACTCA',
                                  'YATGCTGCCTCxCCGTAGGAGT',
                                  'Control',
                                  'ATGACCGATTRGACCAG',
                                  'Control_mouse_I.D._356']]
        expected_comments = [
            'Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).']
        expected_errors = [
            'Treatment found in header 2 times.  Header fields must be unique.\t0,3',
            'Treatment found in header 2 times.  Header fields must be unique.\t0,4',
            'Invalid DNA sequence detected: YATGCTGCCTCxCCGTAGGAGT\t3,2',
            'Invalid DNA sequence detected: AACTCGTCGATGN\t2,1']
        expected_warnings = [
            'Barcode AACTCGTCGATGN differs than length 12\t2,1',
            'Invalid characters found in PC-355\t2,0',
            'Invalid characters found in Cont^^rol\t1,3',
            'Data field outofbounds found after Description column\t2,6']

        self.assertEqual(header, expected_header)
        self.assertEqual(mapping_data, expected_mapping_data)
        self.assertEqual(comments, expected_comments)
        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_data_fields(self):
        """ Overall data fields check returns expected results """

        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s1', 'ACGT', 'AAAA', 's1_data'],
                        ['s2', 'CGTA', 'AAAA', 's2_data']]
        errors = []
        warnings = []

        errors, warnings = check_data_fields(header,
                                             mapping_data, errors, warnings)

        expected_errors = []
        expected_warnings = []

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_data_fields_errors_warnings(self):
        """ Overall data fields check returns expected results """

        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s1', 'ACGT', 'AAAxA', 's1_data'],
                        ['s_2', 'CGTA', 'AAAA', 's2_data']]
        errors = []
        warnings = []

        errors, warnings = check_data_fields(header,
                                             mapping_data, errors, warnings)

        expected_errors = ['Invalid DNA sequence detected: AAAxA\t1,2']
        expected_warnings = ['Invalid characters found in s_2\t2,0']

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_fields_past_bounds(self):
        """ Flags fields found past the Description column as warnings """

        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s1', 'ACGT', 'AAAA', 's1_data'],
                        ['s2', 'CGTA', 'AAAA', 's2_data', 'out_of_bounds']]
        warnings = []

        warnings = check_fields_past_bounds(header, mapping_data, warnings)

        expected_warnings =\
            ['Data field out_of_bounds found after Description column\t2,4']

        self.assertEqual(warnings, expected_warnings)

    def test_check_chars_data_fields(self):
        """ Flags fields with invalid characters as warnings """

        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', 's2_data']]
        warnings = []

        warnings = check_chars_data_fields(header, mapping_data, warnings)

        expected_warnings = ['Invalid characters found in s-1\t1,0',
                             'Invalid characters found in s1&data\t1,3']

        self.assertEqual(warnings, expected_warnings)

    def test_check_dna_chars_primers(self):
        """ Flags primer fields with invalid characters as errors """

        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 's1&data'],
                        ['s2', 'CGTA', 'AAA1A', 's2_data']]
        errors = []

        errors = check_dna_chars_primers(header, mapping_data, errors)

        expected_errors = ['Invalid DNA sequence detected: AAA1A\t2,2']

        self.assertEqual(errors, expected_errors)

        # Should be able to suppress LinkerPrimerSequence check, won't
        # suppress ReversePrimer check

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
             'ReversePrimer', 'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 'ACGT', 's1&data'],
                        ['s2', 'CGTA', 'AAA1A', 'ACGTF', 's2_data']]
        errors = []

        errors = check_dna_chars_primers(header, mapping_data, errors,
                                         disable_primer_check=True)

        expected_errors = ['Invalid DNA sequence detected: ACGTF\t2,3']

        self.assertEqual(errors, expected_errors)

    def test_check_dna_chars_bcs(self):
        """ Flags barcode fields with invalid characters as errors """

        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 's1&data'],
                        ['s2', 'CGTA', 'AAA1A', 's2_data']]
        errors = []

        errors = check_dna_chars_bcs(header, mapping_data, errors)

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should find no errors

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
             'ReversePrimer', 'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 'ACGT', 's1&data'],
                        ['s2', 'C1GTA', 'AAA1A', 'ACGTF', 's2_data']]
        errors = []

        errors = check_dna_chars_bcs(header, mapping_data, errors,
                                     has_barcodes=False)

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should find errors with has_barcodes=True

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
             'ReversePrimer', 'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 'ACGT', 's1&data'],
                        ['s2', 'CNGTA', 'AAA1A', 'ACGTF', 's2_data']]
        errors = []

        errors = check_dna_chars_bcs(header, mapping_data, errors,
                                     has_barcodes=True)

        expected_errors = ['Invalid DNA sequence detected: CNGTA\t2,1']

        self.assertEqual(errors, expected_errors)

    def test_check_bcs_lengths(self):
        """ Flags barcodes of different lengths as warnings """

        # should not give any warnings for equal lengths
        header =\
            ['SampleID',
             'BarcodeSequence',
             'LinkerPrimerSequence',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 's1&data'],
                        ['s2', 'CGTA', 'AAA1A', 's2_data']]
        warnings = []

        warnings = check_bcs_lengths(header, mapping_data, warnings)

        expected_warnings = []

        self.assertEqual(warnings, expected_warnings)

        # Should give warning for different lengths

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
             'ReversePrimer', 'Description']
        mapping_data = [['s-1', 'ACGT', 'AARNCWSVDAA', 'ACGT', 's1&data'],
                        ['s2', 'C1GTA', 'AAA1A', 'ACGTF', 's2_data'],
                        ['s3', 'CGTA', 'AAAAA', 'ACGTF', 's3_data']]
        warnings = []

        warnings = check_bcs_lengths(header, mapping_data, warnings)

        expected_warnings = ['Barcode C1GTA differs than length 4\t2,1']

        self.assertEqual(warnings, expected_warnings)

    def test_check_bc_duplicates_default_correct(self):
        """ Handles duplicate checks of barcodes and added demultiplex data """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field=None)

        expected_errors = []

        self.assertEqual(errors, expected_errors)

    def test_check_bc_duplicates_default_dups(self):
        """ Handles duplicate checks of barcodes and added demultiplex data """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'ACGT', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field=None)

        expected_errors = ['Duplicate barcode ACGT found.\t1,1',
                           'Duplicate barcode ACGT found.\t2,1']

        self.assertEqual(errors, expected_errors)

    def test_check_bc_duplicates_disable_bcs_dups(self):
        """ Handles duplicate checks of no barcodes or added demultiplex data """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'ACGT', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=False,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field=None)

        # Should return an error if more than one SampleID present and
        # barcodes disabled and no added_demultiplex_field
        expected_errors = [
            "If no barcodes are present, and the added_demultiplex_field option isn't used, only a single SampleID can be present.\t-1,-1"]

        self.assertEqual(errors, expected_errors)

        # If a single sample present, should not raise any errors

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=False,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field=None)

        expected_errors = []

        self.assertEqual(errors, expected_errors)

    def test_check_bc_duplicates_var_len_no_dupes(self):
        """ Handles duplicate checks of barcodes and added demultiplex data """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'ACGT', 'TAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=True,
                                     added_demultiplex_field=None)

        # combination of primer seq and barcodes to match largest barcode
        # present is ACGTA and ACGTT, so should not get a duplicate hit.
        expected_errors = []

        self.assertEqual(errors, expected_errors)

    def test_check_bc_duplicates_var_len_dupes(self):
        """ Handles duplicate checks of barcodes and added demultiplex data """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'ACGT', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=True,
                                     added_demultiplex_field=None)

        # Barcode 1 is the largest, with 5 nts, is sequence ACGTA.  When the
        # last base at 5' end of primer is added to barcode 2, there is a
        # duplicate, as this is also ACGTA.
        expected_errors = [
            'Duplicate barcode and primer fragment sequence ACGTA found.\t1,1',
            'Duplicate barcode and primer fragment sequence ACGTA found.\t2,1']

        self.assertEqual(errors, expected_errors)

    def test_check_bc_duplicates_added_demultiplex(self):
        """ Handles duplicate checks of barcodes and added demultiplex data """

        # Should not find any duplicates
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should not find any duplicates with var length turned on.
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=True,
                                     added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should not find errors when only looking at added field
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=False,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

    def test_check_bc_duplicates_added_demultiplex_finds_dups(self):
        """ Handles duplicate checks of barcodes and added demultiplex data """

        # Should find duplicates
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'CGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '1', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field='run_prefix')

        expected_errors = [
            'Duplicate barcode and added demultiplex field CGTA1 found.\t1,1',
            'Duplicate barcode and added demultiplex field CGTA1 found.\t2,1']

        self.assertEqual(errors, expected_errors)

        # Should find duplicates with var length turned on
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'CGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTAA', 'AAAA', '1', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=True,
                                     variable_len_barcodes=True,
                                     added_demultiplex_field='run_prefix')

        expected_errors = [
            'Duplicate barcode and added demultiplex field CGTAA1 found.\t1,1',
            'Duplicate barcode and added demultiplex field CGTAA1 found.\t2,1']

        self.assertEqual(errors, expected_errors)

        # Should find duplicates when just using added fields
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'CGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '1', 's2_data']]
        errors = []

        errors = check_bc_duplicates(header,
                                     mapping_data,
                                     errors,
                                     has_barcodes=False,
                                     variable_len_barcodes=False,
                                     added_demultiplex_field='run_prefix')

        expected_errors = [
            'Duplicate added demultiplex field 1 found.\t1,3',
            'Duplicate added demultiplex field 1 found.\t2,3']

        self.assertEqual(errors, expected_errors)

    def test_check_fixed_len_bcs_dups(self):
        """ Properly detects duplicates of fixed length barcodes """

        # Should not find any duplicates
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_fixed_len_bcs_dups(header,
                                          mapping_data,
                                          errors)

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should find duplicates

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'CGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_fixed_len_bcs_dups(header,
                                          mapping_data,
                                          errors)

        expected_errors = [
            'Duplicate barcode CGTA found.\t1,1',
            'Duplicate barcode CGTA found.\t2,1']

        self.assertEqual(errors, expected_errors)

    def test_check_fixed_len_bcs_dups_mixed_caps(self):

        # Should find duplicates with mixed caps

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'CGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'cgta', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_fixed_len_bcs_dups(header,
                                          mapping_data,
                                          errors)

        expected_errors = [
            'Duplicate barcode CGTA found.\t1,1',
            'Duplicate barcode CGTA found.\t2,1']

        self.assertEqual(errors, expected_errors)

    def test_check_variable_len_bcs_dups(self):
        """ Ensures that slices of barcodes + 5' primer fragments are unique """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'ACGT', 'TAAA', '2', 's2_data']]
        errors = []

        errors = check_variable_len_bcs_dups(header,
                                             mapping_data,
                                             errors)

        # combination of primer seq and barcodes to match largest barcode
        # present is ACGTA and ACGTT, so should not get a duplicate hit.
        expected_errors = []

        self.assertEqual(errors, expected_errors)

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGTA', 'AAAA', '1', 's1&data'],
                        ['s2', 'ACGT', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_variable_len_bcs_dups(header,
                                             mapping_data,
                                             errors)

        # Barcode 1 is the largest, with 5 nts, is sequence ACGTA.  When the
        # last base at 5' end of primer is added to barcode 2, there is a
        # duplicate, as this is also ACGTA.
        expected_errors = [
            'Duplicate barcode and primer fragment sequence ACGTA found.\t1,1',
            'Duplicate barcode and primer fragment sequence ACGTA found.\t2,1']

        self.assertEqual(errors, expected_errors)

    def test_check_added_demultiplex_dups(self):
        """ Checks barcodes/added demultiplex field for unique combinations """

        # Should not find any duplicates
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_added_demultiplex_dups(header, mapping_data, errors,
                                              has_barcodes=True, added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should not find any duplicates with var length turned on.
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_added_demultiplex_dups(header,
                                              mapping_data,
                                              errors,
                                              has_barcodes=True,
                                              added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should not find errors when only looking at added field
        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_added_demultiplex_dups(header,
                                              mapping_data,
                                              errors,
                                              has_barcodes=False,
                                              added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

    def test_check_sampleid_duplicates(self):
        """ Checks that all SampleIDs are unique """

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s2', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_sampleid_duplicates(header, mapping_data, errors)
        # Should not find duplicates
        expected_errors = []

        self.assertEqual(errors, expected_errors)

        header =\
            ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
             'Description']
        mapping_data = [['s-1', 'ACGT', 'AAAA', '1', 's1&data'],
                        ['s-1', 'CGTA', 'AAAA', '2', 's2_data']]
        errors = []

        errors = check_sampleid_duplicates(header, mapping_data, errors)
        # Should find duplicates
        expected_errors = [
            'Duplicate SampleID s-1 found.\t1,0',
            'Duplicate SampleID s-1 found.\t2,0']

        self.assertEqual(errors, expected_errors)

    def test_check_header(self):
        """ All header problems are detected """

        # Default header, should not generate any errors/warnings
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'Description']
        errors = []
        warnings = []

        errors, warnings = check_header(header,
                                        errors,
                                        warnings,
                                        sample_id_ix=0,
                                        desc_ix=3,
                                        bc_ix=1,
                                        linker_primer_ix=2,
                                        added_demultiplex_field=None)

        expected_errors = []
        expected_warnings = []

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_header_missing_fields(self):
        """ All header problems are detected """

        # Default header, should not generate any errors/warnings
        header = ['AAA', 'XXX', 'YYY',
                  'ZZZ']
        errors = []
        warnings = []

        errors, warnings = check_header(header,
                                        errors,
                                        warnings,
                                        sample_id_ix=0,
                                        desc_ix=3,
                                        bc_ix=1,
                                        linker_primer_ix=2,
                                        added_demultiplex_field=None)

        expected_errors = [
            'Found header field AAA, expected field SampleID\t0,0',
            'Found header field XXX, expected field BarcodeSequence\t0,1',
            'Found header field YYY, expected field LinkerPrimerSequence\t0,2',
            'Found header field ZZZ, last field should be Description\t0,3']
        expected_warnings = []

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_header_bad_chars(self):
        """ All header problems are detected """

        # Default header, should not generate any errors/warnings
        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'AAA.^^2',
            'Description']
        errors = []
        warnings = []

        errors, warnings = check_header(header,
                                        errors,
                                        warnings,
                                        sample_id_ix=0,
                                        desc_ix=4,
                                        bc_ix=1,
                                        linker_primer_ix=2,
                                        added_demultiplex_field=None)

        expected_errors = []
        expected_warnings = [
            'Found invalid character in AAA.^^2 header field.\t0,3']

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_header_missing_added_demultiplex(self):
        """ All header problems are detected """

        # Default header, should not generate any errors/warnings
        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
            'Description']
        errors = []
        warnings = []

        errors, warnings = check_header(header,
                                        errors,
                                        warnings,
                                        sample_id_ix=0,
                                        desc_ix=4,
                                        bc_ix=1,
                                        linker_primer_ix=2,
                                        added_demultiplex_field='run_prefix')

        expected_errors = []
        expected_warnings = []

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

        # Default header, should not generate any errors/warnings
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'Description']
        errors = []
        warnings = []

        errors, warnings = check_header(header,
                                        errors,
                                        warnings,
                                        sample_id_ix=0,
                                        desc_ix=3,
                                        bc_ix=1,
                                        linker_primer_ix=2,
                                        added_demultiplex_field='run_prefix')

        expected_errors = ['Missing added demultiplex field run_prefix\t-1,-1']
        expected_warnings = []

        self.assertEqual(errors, expected_errors)
        self.assertEqual(warnings, expected_warnings)

    def test_check_header_dups(self):
        """ Flags header duplicates as errors """

        # Default header, should not generate any errors/warnings
        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
            'Description']
        errors = []

        errors = check_header_dups(header, errors)

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should give errors with dups
        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix', 'run_prefix',
            'Description']
        errors = []

        errors = check_header_dups(header, errors)

        expected_errors = [
            'run_prefix found in header 2 times.  Header fields must be unique.\t0,3',
            'run_prefix found in header 2 times.  Header fields must be unique.\t0,4']

        self.assertEqual(errors, expected_errors)

    def test_check_header_chars(self):
        """ Flags invalid characters in header as warnings """

        # Default header, should not generate any errors/warnings
        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_p-%efix',
            'Description']
        warnings = []

        warnings = check_header_chars(header, warnings)

        expected_warnings = [
            'Found invalid character in run_p-%efix header field.\t0,3']

        self.assertEqual(warnings, expected_warnings)

    def test_check_header_required_fields(self):
        """ Flags missing header fields as errors """

        # Default header, should not generate any errors/warnings
        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
            'Description']
        errors = []

        errors = check_header_required_fields(header,
                                              errors,
                                              sample_id_ix=0,
                                              desc_ix=4,
                                              bc_ix=1,
                                              linker_primer_ix=2,
                                              added_demultiplex_field='run_prefix')

        expected_errors = []

        self.assertEqual(errors, expected_errors)

        # Should find all as errors if not named correctly
        header = ['AAA', 'BBB', 'CCC', 'DDD',
                  'EEE']
        errors = []

        errors = check_header_required_fields(header,
                                              errors,
                                              sample_id_ix=0,
                                              desc_ix=4,
                                              bc_ix=1,
                                              linker_primer_ix=2,
                                              added_demultiplex_field='run_prefix')

        expected_errors = [
            'Found header field AAA, expected field SampleID\t0,0',
            'Found header field BBB, expected field BarcodeSequence\t0,1',
            'Found header field CCC, expected field LinkerPrimerSequence\t0,2',
            'Found header field EEE, last field should be Description\t0,4',
            'Missing added demultiplex field run_prefix\t-1,-1']

        self.assertEqual(errors, expected_errors)

    def test_correct_mapping_data(self):
        """ Properly replaces invalid characters in mapping data """

        header = [
            'SampleID', 'BarcodeSequence', 'LinkerPrimerSequence', 'run_prefix',
            'Description']
        mapping_data = [
            ['AA_2', 'ACCAAGGACTT', 'ACGGATACCGAGx', '26^3', 'A-2'],
            ['AA.3', 'ACAGGATAAC', 'ARRACGGA', '12__', 'AA%3']]

        corrected_data = correct_mapping_data(mapping_data, header)

        # Should replace SampleID invalid characters with periods, other
        # invalid characters with underscores, should not change DNA fields.
        expected_corrected_data = [['AA.2',
                                    'ACCAAGGACTT',
                                    'ACGGATACCGAGx',
                                    '26_3',
                                    'A-2'],
                                   ['AA.3',
                                    'ACAGGATAAC',
                                    'ARRACGGA',
                                    '12__',
                                    'AA%3']]

        self.assertEqual(corrected_data, expected_corrected_data)


# Input data
sample_correct_mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	ReversePrimer	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356"""

sample_errors_mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	ReversePrimer	NotDescription
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.355	AGCACGAGCCxTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._355
 PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356"""


sample_warnings_mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatm-ent	ReversePrimer	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC_355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Co&ntrol	ATGACCGATTRGACCAG	Control_mouse_I.D._355	OutOfBounds
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356"""

sample_errors_warnings_mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	Treatment	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Cont^^rol	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC-355	AACTCGTCGATGN	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._355	outofbounds
PC.356	ACAGACCACTCA	YATGCTGCCTCxCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356"""

empty_fields_mapping_data = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	ReversePrimer	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356"""

# Expected output data
expected_html_data_correct_input = """<html>\n<head>\n\n<script type="text/javascript" src="./overlib.js"></script>\n</head>\n<body bgcolor="white"> <h1>No errors or warnings detected.<br></h1><h1>Mapping file error and warning details.</h1>\nNotes for interpreting this report:\n<ul>\n    <li>Errors will be listed in red, warnings in yellow.\n    <li>Mouse over an error or warning in a cell for more details.\n    <li>Errors in the header row may mask other errors, so these should be corrected first.\n    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>validate_mapping_file.py</tt> until no errors (nor warnings, ideally) are found.\n</ul>\n<p>\nSome general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):\n<ul>\n    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.\n    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>\n    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.\n</ul>\nGeneral issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr></tr></table><br>\n<table border="2" cellspacing="0" cellpadding="5">\n\n<tr></tr>\n<tr>\n<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Treatment</th><th>ReversePrimer</th><th>Description</th>\n</tr>\n\n<tr>\n<tr><th><tt>PC.354</tt></th><th><tt>AGCACGAGCCTA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._354\n</tt></th></tr><tr><th><tt>PC.355</tt></th><th><tt>AACTCGTCGATG</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._355\n</tt></th></tr><tr><th><tt>PC.356</tt></th><th><tt>ACAGACCACTCA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._356</tt></th></tr>\n</tr>\n</table>\n\n</body>\n</html>"""
expected_corrected_data_correct_input = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	ReversePrimer	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356
"""

expected_log_data_correct_input = """No errors or warnings found in mapping file."""

expected_html_errors_output = """<html>\n<head>\n\n<script type="text/javascript" src="./overlib.js"></script>\n</head>\n<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>\nNotes for interpreting this report:\n<ul>\n    <li>Errors will be listed in red, warnings in yellow.\n    <li>Mouse over an error or warning in a cell for more details.\n    <li>Errors in the header row may mask other errors, so these should be corrected first.\n    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>validate_mapping_file.py</tt> until no errors (nor warnings, ideally) are found.\n</ul>\n<p>\nSome general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):\n<ul>\n    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.\n    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>\n    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.\n</ul>\nGeneral issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr></tr></table><br>\n<table border="2" cellspacing="0" cellpadding="5">\n\n<tr></tr>\n<tr>\n<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Treatment</th><th>ReversePrimer</th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Found header field NotDescription, last field should be Description<br>\');" onmouseout="return nd();"><font color=white>NotDescription</a></th>\n</tr>\n\n<tr>\n<tr><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Duplicate SampleID PC.355 found.<br>Location (SampleID,Header Field)<br>PC.355,SampleID\');" onmouseout="return nd();"><font color=white><tt>PC.355</tt></a></th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Barcode AGCACGAGCCxTA differs than length 12<br>Invalid DNA sequence detected: AGCACGAGCCxTA<br>Location (SampleID,Header Field)<br>PC.355,BarcodeSequence\');" onmouseout="return nd();"><font color=white><tt>AGCACGAGCCxTA</tt></a></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._354\n</tt></th></tr><tr><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Duplicate SampleID PC.355 found.<br>Location (SampleID,Header Field)<br>PC.355,SampleID\');" onmouseout="return nd();"><font color=white><tt>PC.355</tt></a></th><th><tt>AACTCGTCGATG</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._355\n</tt></th></tr><tr><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid characters found in  PC.356<br>Location (SampleID,Header Field)<br> PC.356,SampleID\');" onmouseout="return nd();"><font color=black><tt> PC.356</tt></a></th><th><tt>ACAGACCACTCA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._356</tt></th></tr>\n</tr>\n</table>\n\n</body>\n</html>"""
expected_data_errors_corrected_output = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tReversePrimer\tNotDescription\n#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).\nPC.355\tAGCACGAGCCxTA\tYATGCTGCCTCCCGTAGGAGT\tControl\tATGACCGATTRGACCAG\tControl_mouse_I.D._354\nPC.355\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\tATGACCGATTRGACCAG\tControl_mouse_I.D._355\n.PC.356\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tControl\tATGACCGATTRGACCAG\tControl_mouse_I.D._356\n"""

expected_data_log_errors_output = """# Errors and warnings are written as a tab separated columns, with the first column showing the error or warning, and the second column contains the location of the error or warning, written as row,column, where 0,0 is the top left header item (SampleID).  Problems not specific to a particular data cell will be listed as having 'no location'.\nErrors -----------------------------\nFound header field NotDescription, last field should be Description\t0,5\nInvalid DNA sequence detected: AGCACGAGCCxTA\t1,1\nDuplicate SampleID PC.355 found.\t1,0\nDuplicate SampleID PC.355 found.\t2,0\nWarnings ---------------------------\nBarcode AGCACGAGCCxTA differs than length 12\t1,1\nInvalid characters found in  PC.356\t3,0\n"""

expected_html_errors_suppressed_bcs = """<html>\n<head>\n\n<script type="text/javascript" src="./overlib.js"></script>\n</head>\n<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>\nNotes for interpreting this report:\n<ul>\n    <li>Errors will be listed in red, warnings in yellow.\n    <li>Mouse over an error or warning in a cell for more details.\n    <li>Errors in the header row may mask other errors, so these should be corrected first.\n    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>validate_mapping_file.py</tt> until no errors (nor warnings, ideally) are found.\n</ul>\n<p>\nSome general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):\n<ul>\n    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.\n    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>\n    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.\n</ul>\nGeneral issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr><td bgcolor="red"><font color="white">If no barcodes are present, and the added_demultiplex_field option isn\'t used, only a single SampleID can be present.<font color="black"></td></tr></table><br>\n<table border="2" cellspacing="0" cellpadding="5">\n\n<tr></tr>\n<tr>\n<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th>Treatment</th><th>ReversePrimer</th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Found header field NotDescription, last field should be Description<br>\');" onmouseout="return nd();"><font color=white>NotDescription</a></th>\n</tr>\n\n<tr>\n<tr><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Duplicate SampleID PC.355 found.<br>Location (SampleID,Header Field)<br>PC.355,SampleID\');" onmouseout="return nd();"><font color=white><tt>PC.355</tt></a></th><th><tt>AGCACGAGCCxTA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._354\n</tt></th></tr><tr><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Duplicate SampleID PC.355 found.<br>Location (SampleID,Header Field)<br>PC.355,SampleID\');" onmouseout="return nd();"><font color=white><tt>PC.355</tt></a></th><th><tt>AACTCGTCGATG</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._355\n</tt></th></tr><tr><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid characters found in  PC.356<br>Location (SampleID,Header Field)<br> PC.356,SampleID\');" onmouseout="return nd();"><font color=black><tt> PC.356</tt></a></th><th><tt>ACAGACCACTCA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._356</tt></th></tr>\n</tr>\n</table>\n\n</body>\n</html>"""
expected_output_log_errors_bcs_suppressed = """# Errors and warnings are written as a tab separated columns, with the first column showing the error or warning, and the second column contains the location of the error or warning, written as row,column, where 0,0 is the top left header item (SampleID).  Problems not specific to a particular data cell will be listed as having 'no location'.\nErrors -----------------------------\nFound header field NotDescription, last field should be Description\t0,5\nIf no barcodes are present, and the added_demultiplex_field option isn't used, only a single SampleID can be present.\tno location\nDuplicate SampleID PC.355 found.\t1,0\nDuplicate SampleID PC.355 found.\t2,0\nWarnings ---------------------------\nInvalid characters found in  PC.356\t3,0\n"""

expected_html_output_warnings = """<html>\n<head>\n\n<script type="text/javascript" src="./overlib.js"></script>\n</head>\n<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>\nNotes for interpreting this report:\n<ul>\n    <li>Errors will be listed in red, warnings in yellow.\n    <li>Mouse over an error or warning in a cell for more details.\n    <li>Errors in the header row may mask other errors, so these should be corrected first.\n    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>validate_mapping_file.py</tt> until no errors (nor warnings, ideally) are found.\n</ul>\n<p>\nSome general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):\n<ul>\n    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.\n    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>\n    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.\n</ul>\nGeneral issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr></tr></table><br>\n<table border="2" cellspacing="0" cellpadding="5">\n\n<tr></tr>\n<tr>\n<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Found invalid character in Treatm-ent header field.<br>\');" onmouseout="return nd();"><font color=black>Treatm-ent</a></th><th>ReversePrimer</th><th>Description</th>\n</tr>\n\n<tr>\n<tr><th><tt>PC.354</tt></th><th><tt>AGCACGAGCCTA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._354\n</tt></th></tr><tr><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid characters found in PC_355<br>Location (SampleID,Header Field)<br>PC_355,SampleID\');" onmouseout="return nd();"><font color=black><tt>PC_355</tt></a></th><th><tt>AACTCGTCGATG</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid characters found in Co&ntrol<br>Location (SampleID,Header Field)<br>PC_355,Treatm-ent\');" onmouseout="return nd();"><font color=black><tt>Co&ntrol</tt></a></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._355</tt></th><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Data field OutOfBounds found after Description column<br>Location (SampleID,Header Field)<br>PC_355,no header\');" onmouseout="return nd();"><font color=black><tt>OutOfBounds\n</tt></a></th></tr><tr><th><tt>PC.356</tt></th><th><tt>ACAGACCACTCA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._356</tt></th></tr>\n</tr>\n</table>\n\n</body>\n</html>"""

expected_corrected_warnings_output = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatm-ent	ReversePrimer	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Co_ntrol	ATGACCGATTRGACCAG	Control_mouse_I.D._355	OutOfBounds
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356
"""

expected_log_warnings_output = """# Errors and warnings are written as a tab separated columns, with the first column showing the error or warning, and the second column contains the location of the error or warning, written as row,column, where 0,0 is the top left header item (SampleID).  Problems not specific to a particular data cell will be listed as having 'no location'.
Errors -----------------------------
Warnings ---------------------------
Found invalid character in Treatm-ent header field.	0,3
Invalid characters found in PC_355	2,0
Invalid characters found in Co&ntrol	2,3
Data field OutOfBounds found after Description column	2,6
"""

expected_html_errors_warnings_output = """<html>\n<head>\n\n<script type="text/javascript" src="./overlib.js"></script>\n</head>\n<body bgcolor="white"> <h1>Mapping file error and warning details.</h1>\nNotes for interpreting this report:\n<ul>\n    <li>Errors will be listed in red, warnings in yellow.\n    <li>Mouse over an error or warning in a cell for more details.\n    <li>Errors in the header row may mask other errors, so these should be corrected first.\n    <li>Modifications to your mapping file to fix certain issues may result in different errors. You should run <tt>validate_mapping_file.py</tt> until no errors (nor warnings, ideally) are found.\n</ul>\n<p>\nSome general rules about formatting mapping files (see <a href="http://qiime.org/documentation/file_formats.html#metadata-mapping-files">here</a> for additional details):\n<ul>\n    <li>Header characters should only contain alphanumeric and <tt>_</tt> characters only.\n    <li>Valid characters for SampleID fields are alphanumeric and <tt>.</tt> only.<br>\n    <li>Other fields allow alphanumeric and <tt>+-%./ :,;_</tt> characters.\n</ul>\nGeneral issues with your mapping file (i.e., those that do not pertain to a particular cell) will be listed here, if any:<table border="1" cellspacing="0" cellpadding="7"><tr><td bgcolor="red"><font color="white">Missing added demultiplex field DoesNotExist<font color="black"></td></tr></table><br>\n<table border="2" cellspacing="0" cellpadding="5">\n\n<tr></tr>\n<tr>\n<th>SampleID</th><th>BarcodeSequence</th><th>LinkerPrimerSequence</th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Treatment found in header 2 times.  Header fields must be unique.<br>\');" onmouseout="return nd();"><font color=white>Treatment</a></th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Treatment found in header 2 times.  Header fields must be unique.<br>\');" onmouseout="return nd();"><font color=white>Treatment</a></th><th>Description</th>\n</tr>\n\n<tr>\n<tr><th><tt>PC.354</tt></th><th><tt>AGCACGAGCCTA</tt></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid characters found in Cont^^rol<br>Location (SampleID,Header Field)<br>PC.354,Treatment\');" onmouseout="return nd();"><font color=black><tt>Cont^^rol</tt></a></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._354\n</tt></th></tr><tr><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid characters found in PC-355<br>Location (SampleID,Header Field)<br>PC-355,SampleID\');" onmouseout="return nd();"><font color=black><tt>PC-355</tt></a></th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Barcode AACTCGTCGATGN differs than length 12<br>Invalid DNA sequence detected: AACTCGTCGATGN<br>Location (SampleID,Header Field)<br>PC-355,BarcodeSequence\');" onmouseout="return nd();"><font color=white><tt>AACTCGTCGATGN</tt></a></th><th><tt>YATGCTGCCTCCCGTAGGAGT</tt></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._355</tt></th><th bgcolor=yellow><a href="javascript:void(0);" onmouseover="return overlib(\'Data field outofbounds found after Description column<br>Location (SampleID,Header Field)<br>PC-355,no header\');" onmouseout="return nd();"><font color=black><tt>outofbounds\n</tt></a></th></tr><tr><th><tt>PC.356</tt></th><th><tt>ACAGACCACTCA</tt></th><th bgcolor=red><a href="javascript:void(0);" onmouseover="return overlib(\'Invalid DNA sequence detected: YATGCTGCCTCxCCGTAGGAGT<br>Location (SampleID,Header Field)<br>PC.356,LinkerPrimerSequence\');" onmouseout="return nd();"><font color=white><tt>YATGCTGCCTCxCCGTAGGAGT</tt></a></th><th><tt>Control</tt></th><th><tt>ATGACCGATTRGACCAG</tt></th><th><tt>Control_mouse_I.D._356</tt></th></tr>\n</tr>\n</table>\n\n</body>\n</html>"""

expected_corrected_data_errors_warnings = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	Treatment	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Cont__rol	ATGACCGATTRGACCAG	Control_mouse_I.D._354
PC.355	AACTCGTCGATGN	YATGCTGCCTCCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._355	outofbounds
PC.356	ACAGACCACTCA	YATGCTGCCTCxCCGTAGGAGT	Control	ATGACCGATTRGACCAG	Control_mouse_I.D._356
"""

expected_log_errors_warnings_output = """# Errors and warnings are written as a tab separated columns, with the first column showing the error or warning, and the second column contains the location of the error or warning, written as row,column, where 0,0 is the top left header item (SampleID).  Problems not specific to a particular data cell will be listed as having 'no location'.
Errors -----------------------------
Treatment found in header 2 times.  Header fields must be unique.	0,3
Treatment found in header 2 times.  Header fields must be unique.	0,4
Missing added demultiplex field DoesNotExist	no location
Invalid DNA sequence detected: YATGCTGCCTCxCCGTAGGAGT	3,2
Invalid DNA sequence detected: AACTCGTCGATGN	2,1
Warnings ---------------------------
Barcode AACTCGTCGATGN differs than length 12	2,1
Invalid characters found in PC-355	2,0
Invalid characters found in Cont^^rol	1,3
Data field outofbounds found after Description column	2,6
"""

if __name__ == '__main__':
    main()
