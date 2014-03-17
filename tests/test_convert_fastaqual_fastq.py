#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka, Abhisaar Yadav"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka", "Abhisaar Yadav", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

# Reviewed by William Walters

from os.path import sep, split, splitext, exists, join
from shutil import rmtree
from os import chmod

from unittest import TestCase, main
from cogent.util.misc import remove_files, get_random_directory_name

from qiime.util import get_tmp_filename, create_dir

from qiime.convert_fastaqual_fastq import (convert_fastq, convert_fastaqual,
                                           convert_fastaqual_fastq,
                                           get_filename_with_new_ext)


class MakeFastqTests(TestCase):

    """ Unit tests for the convert_fastaqual_fastq.py module """

    def setUp(self):
        self._files_to_remove = []

        self.qual_file_path = get_tmp_filename(prefix='qual_', suffix='.qual')
        self.fasta_file_path = get_tmp_filename(prefix='fasta_', suffix='.fna')
        self.nolabel_qual_file_path = get_tmp_filename(prefix='qual_',
                                                       suffix='.qual')
        self.noseq_qual_file_path = get_tmp_filename(prefix='qual_',
                                                     suffix='.qual')

        qual_file = open(self.qual_file_path, 'w')
        fasta_file = open(self.fasta_file_path, 'w')

        qual_file.write(qual_test_string)
        qual_file.close()

        fasta_file.write(fasta_test_string)
        fasta_file.close()

        # Error testing files
        nolabel_qual_file = open(self.nolabel_qual_file_path, 'w')
        nolabel_qual_file.write(nolabel_qual_test_string)
        nolabel_qual_file.close()
        noseq_qual_file = open(self.noseq_qual_file_path, 'w')
        noseq_qual_file.write(noseq_qual_test_string)
        noseq_qual_file.close()
        self.read_only_output_dir = get_tmp_filename(prefix='read_only_',
                                                     suffix='/')
        create_dir(self.read_only_output_dir)
        # Need read only directory to test errors for files written during
        # fastq/fasta iteration.
        chmod(self.read_only_output_dir, 0o555)

        self.output_dir = get_tmp_filename(prefix='convert_fastaqual_fastq_',
                                           suffix='/')
        self.output_dir += sep

        create_dir(self.output_dir)

        self._files_to_remove.append(self.qual_file_path)
        self._files_to_remove.append(self.fasta_file_path)

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
        if exists(self.read_only_output_dir):
            rmtree(self.read_only_output_dir)

    def test_default_settings(self):
        """ Handles conversions with default settings """
        convert_fastq(self.fasta_file_path, self.qual_file_path,
                      output_directory=self.output_dir)

        actual_output_file_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.fastq',
            self.output_dir)

        actual_output_file = open(actual_output_file_path)
        actual_output = actual_output_file.read()
        actual_output_file.close()
        self._files_to_remove.append(actual_output_file_path)

        self.assertEquals(actual_output, expected_fastq_default_options)

    def test_full_fasta_headers(self):
        """ Properly retains full fasta headers """
        convert_fastq(self.fasta_file_path, self.qual_file_path,
                      full_fasta_headers=True, output_directory=self.output_dir)

        actual_output_file_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.fastq',
            self.output_dir)

        actual_output_file = open(actual_output_file_path)
        actual_output = actual_output_file.read()
        actual_output_file.close()
        self._files_to_remove.append(actual_output_file_path)

        self.assertEquals(actual_output, expected_fastq_full_fasta_headers)

    def test_full_fasta_full_fastq(self):
        """ Properly writes labels to quality score headers """
        convert_fastq(self.fasta_file_path, self.qual_file_path,
                      full_fasta_headers=True, full_fastq=True,
                      output_directory=self.output_dir)

        actual_output_file_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.fastq',
            self.output_dir)

        actual_output_file = open(actual_output_file_path)
        actual_output = actual_output_file.read()
        actual_output_file.close()
        self._files_to_remove.append(actual_output_file_path)

        self.assertEquals(actual_output, expected_fastq_full_fasta_full_fastq)

    def test_multiple_output_files(self):
        """ properly writes multiple fasta files for each sampleID"""
        convert_fastq(self.fasta_file_path, self.qual_file_path,
                      multiple_output_files=True,
                      output_directory=self.output_dir,
                      per_file_buffer_size=23)

        sample_ids = [('PC.634', expected_fastq_634_default),
                      ('PC.354', expected_fastq_354_default),
                      ('PC.481', expected_fastq_481_default)]
        for sample_id, expected_output in sample_ids:
            actual_output_file_path = get_filename_with_new_ext(
                self.fasta_file_path,
                '_' + sample_id + '.fastq',
                self.output_dir)

            actual_output_file = open(actual_output_file_path)
            actual_output = actual_output_file.read()
            actual_output_file.close()
            self._files_to_remove.append(actual_output_file_path)

            self.assertEquals(actual_output, expected_output)

    def test_ascii_increment(self):
        """ Tests for proper range of ascii increments """
        self.assertRaises(ValueError, convert_fastq, self.fasta_file_path,
                          self.qual_file_path, ascii_increment=140, output_directory=self.output_dir)
        self.assertRaises(ValueError, convert_fastq, self.fasta_file_path,
                          self.qual_file_path, ascii_increment=10, output_directory=
                          self.output_dir)

    def test_fastq_output(self):
        """ Raises errors when can't write output file """
        self.assertRaises(IOError, convert_fastq, self.fasta_file_path,
                          self.qual_file_path, output_directory=self.read_only_output_dir)

    def test_qual_label(self):
        """ Raises error if mismatch between quality and fasta labels """
        self.assertRaises(KeyError, convert_fastq, self.fasta_file_path,
                          self.nolabel_qual_file_path, output_directory=self.output_dir)

    def test_qual_seq_length(self):
        """ Raises error if mismatch between length of fasta and qual scores """
        self.assertRaises(KeyError, convert_fastq, self.fasta_file_path,
                          self.noseq_qual_file_path, output_directory=self.output_dir)


class MakeFastaqualTests(TestCase):

    """ Unit tests for the convert_fastaqual_fastq.py module """

    def setUp(self):
        self._files_to_remove = []

        self.fasta_file_path = get_tmp_filename(prefix='fastq_',
                                                suffix='.fastq')

        fastq_file = open(self.fasta_file_path, 'w')

        fastq_file.write(fastq_test_string)
        fastq_file.close()

        # Error testing files
        false_fasta_file = '/'
        false_qual_file = '/'
        self.read_only_output_dir = get_tmp_filename(prefix='read_only_',
                                                     suffix='/')
        create_dir(self.read_only_output_dir)
        chmod(self.read_only_output_dir, 0o555)

        self.output_dir = get_tmp_filename(prefix='convert_fastaqual_fastq_',
                                           suffix='/')
        self.output_dir += sep

        create_dir(self.output_dir)

        self._files_to_remove.append(self.fasta_file_path)

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
        if exists(self.read_only_output_dir):
            rmtree(self.read_only_output_dir)

    def test_default_settings(self):
        """ Converting to fasta/qual files handles default settings """
        convert_fastaqual(self.fasta_file_path,
                          output_directory=self.output_dir)

        actual_output_fasta_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.fna',
            self.output_dir)

        actual_output_qual_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.qual',
            self.output_dir)

        actual_output_fasta = open(actual_output_fasta_path)
        actual_output_qual = open(actual_output_qual_path)
        actual_fasta = actual_output_fasta.read()
        actual_output_fasta.close()
        actual_qual = actual_output_qual.read()
        actual_output_qual.close()
        self._files_to_remove.append(actual_output_fasta_path)
        self._files_to_remove.append(actual_output_qual_path)

        self.assertEquals(actual_fasta, expected_fasta_default_options)
        self.assertEquals(actual_qual, expected_qual_default_options)

    def test_full_fasta_headers(self):
        """ Full headers written to fasta/qual files """
        convert_fastaqual(self.fasta_file_path, full_fasta_headers=True,
                          output_directory=self.output_dir)

        actual_output_fasta_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.fna',
            self.output_dir)

        actual_output_qual_path = get_filename_with_new_ext(
            self.fasta_file_path,
            '.qual',
            self.output_dir)

        actual_output_fasta = open(actual_output_fasta_path)
        actual_output_qual = open(actual_output_qual_path)
        actual_fasta = actual_output_fasta.read()
        actual_output_fasta.close()
        actual_qual = actual_output_qual.read()
        actual_output_qual.close()
        self._files_to_remove.append(actual_output_fasta_path)
        self._files_to_remove.append(actual_output_qual_path)

        self.assertEquals(actual_fasta, expected_fasta_full_fasta_headers)
        self.assertEquals(actual_qual, expected_qual_full_fasta_headers)

    def test_multiple_output_files(self):
        """ Creates one file per sampleID for fasta/qual output """
        convert_fastaqual(self.fasta_file_path,
                          multiple_output_files=True,
                          output_directory=self.output_dir,
                          per_file_buffer_size=23)

        sample_id_s = [('PC.634', expected_fasta_634_default,
                        expected_qual_634_default),
                       ('PC.354', expected_fasta_354_default,
                        expected_qual_354_default),
                       ('PC.481', expected_fasta_481_default,
                        expected_qual_481_default)]
        for sample_id, expected_fasta, expected_qual in sample_id_s:
            actual_output_fasta_path = get_filename_with_new_ext(
                self.fasta_file_path,
                '_' + sample_id + '.fna',
                self.output_dir)

            actual_output_qual_path = get_filename_with_new_ext(
                self.fasta_file_path,
                '_' + sample_id + '.qual',
                self.output_dir)

            actual_output_fasta = open(actual_output_fasta_path)
            actual_output_qual = open(actual_output_qual_path)
            actual_fasta = actual_output_fasta.read()
            actual_output_fasta.close()
            actual_qual = actual_output_qual.read()
            actual_output_qual.close()
            self._files_to_remove.append(actual_output_fasta_path)
            self._files_to_remove.append(actual_output_qual_path)

            self.assertEquals(actual_fasta, expected_fasta)
            self.assertEquals(actual_qual, expected_qual)

    def test_ascii_increment(self):
        """ Detects proper range of ascii increment """
        self.assertRaises(ValueError, convert_fastaqual, self.fasta_file_path,
                          ascii_increment=140, output_directory=self.output_dir)

    def test_fastaqual_output(self):
        """ Raises error if cannot open output filepath """
        self.assertRaises(IOError, convert_fastaqual, self.fasta_file_path,
                          output_directory=self.read_only_output_dir)


class ConvertFastaqualTests(TestCase):

    """ Main function for testing input files, calling proper conversion """

    def setUp(self):
        self._files_to_remove = []

        self.qual_file_path = get_tmp_filename(prefix='qual_', suffix='.qual')
        self.fasta_file_path = get_tmp_filename(prefix='fasta_', suffix='.fna')

        qual_file = open(self.qual_file_path, 'w')
        fasta_file = open(self.fasta_file_path, 'w')
        self.read_only_output_dir = get_tmp_filename(prefix='read_only_',
                                                     suffix='/')
        create_dir(self.read_only_output_dir)
        chmod(self.read_only_output_dir, 0o555)

        self.output_dir = get_tmp_filename(prefix='convert_fastaqual_fastq_',
                                           suffix='/')
        self.output_dir += sep

        create_dir(self.output_dir)

        self._files_to_remove.append(self.qual_file_path)
        self._files_to_remove.append(self.fasta_file_path)

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if exists(self.output_dir):
            rmtree(self.output_dir)
        if exists(self.read_only_output_dir):
            rmtree(self.read_only_output_dir)

    def test_fasta_file(self):
        """ Raises error if cannot open fasta file """
        self.assertRaises(IOError, convert_fastaqual_fastq,
                          self.read_only_output_dir, self.qual_file_path)

    def test_qual_file(self):
        """ Raises error if cannot open qual file """
        self.assertRaises(IOError, convert_fastaqual_fastq,
                          self.fasta_file_path, self.read_only_output_dir)

    def test_conversion_type(self):
        """ Raises error if incorrect conversion type used """
        self.assertRaises(ValueError, convert_fastaqual_fastq,
                          self.fasta_file_path, self.qual_file_path, conversion_type='soijdfl',
                          output_directory=self.output_dir)

    def test_get_filename_with_new_ext(self):
        """ Tests proper function of the utility function. """
        test_paths = [('/from/root/test.xxx', 'test.yyy'),
                      ('../relative/path/test.xxx', 'test.yyy'),
                      ('/double/extension/in/filename/test.zzz.xxx',
                       'test.zzz.yyy')]

        for input, exp_output in test_paths:
            exp_output = join(self.output_dir, exp_output)

            self.assertEquals(
                get_filename_with_new_ext(input, '.yyy', self.output_dir),
                exp_output)

fasta_test_string = '''>PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
'''

qual_test_string = '''>PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
40 40 39 39 38 40 40 40 40 40 40 40 37 37 37 37 37 35 35 35 37 37 37 37 37 35 35 35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 35 32 32 32 32 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37 37 37 36 36 31 31 32 32 36 36 36 36 36 36 36 36 36 36 36 28 27 27 27 26 26 26 30 29 30 29 24 24 24 21 15 15 13 13
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37 37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31 31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25 25 25 28
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
34 34 36 37 36 37 37 37 37 37 37 37 37 37 37 37 36 28 28 28 36 36 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 32 32 32 35 35 37 35 32 32 32 37 37 37 37 37 37 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 35 35 35 37 37 37 37 37 37 37 37 37 37 36 36 36 37 36 35 35 35 37 28 28 28 32 35 37 37 37 36 36 36 37 37 37 37 37 37 35 35 35 35 35 37 37 37 37 36 36 36 37 28 28 28 28 35 36 37 37 37 37 37 37 37 37 37 37 36 33 33 32 31 36 36 33 33 27 27 27 36 31 25 25 25 32 36 36 36 36 36 36 36 36 36 36
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27
'''

nolabel_qual_test_string = '''>PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
40 40 39 39 38 40 40 40 40 40 40 40 37 37 37 37 37 35 35 35 37 37 37 37 37 35 35 35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 35 32 32 32 32 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37 37 37 36 36 31 31 32 32 36 36 36 36 36 36 36 36 36 36 36 28 27 27 27 26 26 26 30 29 30 29 24 24 24 21 15 15 13 13
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37 37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31 31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25 25 25 28
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
34 34 36 37 36 37 37 37 37 37 37 37 37 37 37 37 36 28 28 28 36 36 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 32 32 32 35 35 37 35 32 32 32 37 37 37 37 37 37 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 35 35 35 37 37 37 37 37 37 37 37 37 37 36 36 36 37 36 35 35 35 37 28 28 28 32 35 37 37 37 36 36 36 37 37 37 37 37 37 35 35 35 35 35 37 37 37 37 36 36 36 37 28 28 28 28 35 36 37 37 37 37 37 37 37 37 37 37 36 33 33 32 31 36 36 33 33 27 27 27 36 31 25 25 25 32 36 36 36 36 36 36 36 36 36 36
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27
'''

noseq_qual_test_string = '''>PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
40 40 39 39 38 40 40 40 40 40 40 40 37 37 37 37 37 35 35 35 37 37 37 37 37 35 35 35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 35 32 32 32 32 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37 37 37 36 36 31 31 32 32 36 36 36 36 36 36 36 36 36 36 36 28 27 27 27 26 26 26 30 29 30 29 24 24 24 21 15 15 13 13
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37 37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31 31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25 25 25
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
34 34 36 37 36 37 37 37 37 37 37 37 37 37 37 37 36 28 28 28 36 36 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 32 32 32 35 35 37 35 32 32 32 37 37 37 37 37 37 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 35 35 35 37 37 37 37 37 37 37 37 37 37 36 36 36 37 36 35 35 35 37 28 28 28 32 35 37 37 37 36 36 36 37 37 37 37 37 37 35 35 35 35 35 37 37 37 37 36 36 36 37 28 28 28 28 35 36 37 37 37 37 37 37 37 37 37 37 36 33 33 32 31 36 36 33 33 27 27 27 36 31 25 25 25 32 36 36 36 36 36 36 36 36 36 36
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27
'''

expected_fastq_default_options = '''@PC.634_1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
+
IIHHGIIIIIIIFFFFFDDDFFFFFDDD@@888@666DDFFFEEEEEEFFFFFFFFFFFFFFFFFFFF===EEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFDAAAADFFFFFFFFFFFFFFFFFFEAAAEFDAAAAAAAAEFFFFEE@@AAEEEEEEEEEEE=<<<;;;?>?>999600..
@PC.634_2
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
+
9>BBHHHIHHHIFFFFFFFFFFFAA55555DDFFFFFFFEEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEEEEFFFFFEEEEFFFFFFFFFFFFFFFFEB====EFFFFFFFFFFFFFFFEBBBFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFFFFF===F===FFFFFEEEEE=;;;;=EEEEEEEFGGGGGFFFFF@@@@@@@@@@@@?777::@@@@@@@:::::=
@PC.354_3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
+
CCEFEFFFFFFFFFFFE===EEFFDDDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAAADDFDAAAFFFFFFEEEEEEEEFFFFFFFDDDFFFFFFFFFFFFFFEDDDFFFFFFFFFFEEEFEDDDF===ADFFFEEEFFFFFFDDDDDFFFFEEEF====DEFFFFFFFFFFEBBA@EEBB<<<E@:::AEEEEEEEEEE
@PC.481_4
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
+
IIIIIIIIIIIIFFFFFBBBEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBFFFFFFFFFFFFFFEE====BBBEEFFFFFFFFFFFFFFFFFFFFFFFFEEEEEE@@<<===<6@@EEEEEEEEEEE@@@@@@@
@PC.634_5
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
+
26555@@BDDDDDBBBD@@@DDDDDDD@@@BDDDDDDDDDD@@@;;;;DDDDDDDB@@@DDDDDDDDDDDDDD@@@DDDBBBBBDDDDDDDDDDDDDDDDD?;;;?BDDDDDDDDBBBDB<<:::<
'''

expected_fastq_full_fasta_headers = '''@PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
+
IIHHGIIIIIIIFFFFFDDDFFFFFDDD@@888@666DDFFFEEEEEEFFFFFFFFFFFFFFFFFFFF===EEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFDAAAADFFFFFFFFFFFFFFFFFFEAAAEFDAAAAAAAAEFFFFEE@@AAEEEEEEEEEEE=<<<;;;?>?>999600..
@PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
+
9>BBHHHIHHHIFFFFFFFFFFFAA55555DDFFFFFFFEEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEEEEFFFFFEEEEFFFFFFFFFFFFFFFFEB====EFFFFFFFFFFFFFFFEBBBFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFFFFF===F===FFFFFEEEEE=;;;;=EEEEEEEFGGGGGFFFFF@@@@@@@@@@@@?777::@@@@@@@:::::=
@PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
+
CCEFEFFFFFFFFFFFE===EEFFDDDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAAADDFDAAAFFFFFFEEEEEEEEFFFFFFFDDDFFFFFFFFFFFFFFEDDDFFFFFFFFFFEEEFEDDDF===ADFFFEEEFFFFFFDDDDDFFFFEEEF====DEFFFFFFFFFFEBBA@EEBB<<<E@:::AEEEEEEEEEE
@PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
+
IIIIIIIIIIIIFFFFFBBBEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBFFFFFFFFFFFFFFEE====BBBEEFFFFFFFFFFFFFFFFFFFFFFFFEEEEEE@@<<===<6@@EEEEEEEEEEE@@@@@@@
@PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
+
26555@@BDDDDDBBBD@@@DDDDDDD@@@BDDDDDDDDDD@@@;;;;DDDDDDDB@@@DDDDDDDDDDDDDD@@@DDDBBBBBDDDDDDDDDDDDDDDDD?;;;?BDDDDDDDDBBBDB<<:::<
'''

expected_fastq_full_fasta_full_fastq = '''@PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
+PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
IIHHGIIIIIIIFFFFFDDDFFFFFDDD@@888@666DDFFFEEEEEEFFFFFFFFFFFFFFFFFFFF===EEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFDAAAADFFFFFFFFFFFFFFFFFFEAAAEFDAAAAAAAAEFFFFEE@@AAEEEEEEEEEEE=<<<;;;?>?>999600..
@PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
+PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
9>BBHHHIHHHIFFFFFFFFFFFAA55555DDFFFFFFFEEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEEEEFFFFFEEEEFFFFFFFFFFFFFFFFEB====EFFFFFFFFFFFFFFFEBBBFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFFFFF===F===FFFFFEEEEE=;;;;=EEEEEEEFGGGGGFFFFF@@@@@@@@@@@@?777::@@@@@@@:::::=
@PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
+PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
CCEFEFFFFFFFFFFFE===EEFFDDDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAAADDFDAAAFFFFFFEEEEEEEEFFFFFFFDDDFFFFFFFFFFFFFFEDDDFFFFFFFFFFEEEFEDDDF===ADFFFEEEFFFFFFDDDDDFFFFEEEF====DEFFFFFFFFFFEBBA@EEBB<<<E@:::AEEEEEEEEEE
@PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
+PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
IIIIIIIIIIIIFFFFFBBBEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBFFFFFFFFFFFFFFEE====BBBEEFFFFFFFFFFFFFFFFFFFFFFFFEEEEEE@@<<===<6@@EEEEEEEEEEE@@@@@@@
@PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
+PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
26555@@BDDDDDBBBD@@@DDDDDDD@@@BDDDDDDDDDD@@@;;;;DDDDDDDB@@@DDDDDDDDDDDDDD@@@DDDBBBBBDDDDDDDDDDDDDDDDD?;;;?BDDDDDDDDBBBDB<<:::<
'''

expected_fastq_354_default = '''@PC.354_3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
+
CCEFEFFFFFFFFFFFE===EEFFDDDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAAADDFDAAAFFFFFFEEEEEEEEFFFFFFFDDDFFFFFFFFFFFFFFEDDDFFFFFFFFFFEEEFEDDDF===ADFFFEEEFFFFFFDDDDDFFFFEEEF====DEFFFFFFFFFFEBBA@EEBB<<<E@:::AEEEEEEEEEE
'''

expected_fastq_481_default = '''@PC.481_4
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
+
IIIIIIIIIIIIFFFFFBBBEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBFFFFFFFFFFFFFFEE====BBBEEFFFFFFFFFFFFFFFFFFFFFFFFEEEEEE@@<<===<6@@EEEEEEEEEEE@@@@@@@
'''

expected_fastq_634_default = '''@PC.634_1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
+
IIHHGIIIIIIIFFFFFDDDFFFFFDDD@@888@666DDFFFEEEEEEFFFFFFFFFFFFFFFFFFFF===EEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFDAAAADFFFFFFFFFFFFFFFFFFEAAAEFDAAAAAAAAEFFFFEE@@AAEEEEEEEEEEE=<<<;;;?>?>999600..
@PC.634_2
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
+
9>BBHHHIHHHIFFFFFFFFFFFAA55555DDFFFFFFFEEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEEEEFFFFFEEEEFFFFFFFFFFFFFFFFEB====EFFFFFFFFFFFFFFFEBBBFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFFFFF===F===FFFFFEEEEE=;;;;=EEEEEEEFGGGGGFFFFF@@@@@@@@@@@@?777::@@@@@@@:::::=
@PC.634_5
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
+
26555@@BDDDDDBBBD@@@DDDDDDD@@@BDDDDDDDDDD@@@;;;;DDDDDDDB@@@DDDDDDDDDDDDDD@@@DDDBBBBBDDDDDDDDDDDDDDDDD?;;;?BDDDDDDDDBBBDB<<:::<
'''

fastq_test_string = '''@PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
+
IIHHGIIIIIIIFFFFFDDDFFFFFDDD@@888@666DDFFFEEEEEEFFFFFFFFFFFFFFFFFFFF===EEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFDAAAADFFFFFFFFFFFFFFFFFFEAAAEFDAAAAAAAAEFFFFEE@@AAEEEEEEEEEEE=<<<;;;?>?>999600..
@PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
+
9>BBHHHIHHHIFFFFFFFFFFFAA55555DDFFFFFFFEEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFEEEEEEFFFFFEEEEFFFFFFFFFFFFFFFFEB====EFFFFFFFFFFFFFFFEBBBFFFFFFFFFFFFFFFFFFFFFEEEFFFFFFFFFFFFFFFFFFFF===F===FFFFFEEEEE=;;;;=EEEEEEEFGGGGGFFFFF@@@@@@@@@@@@?777::@@@@@@@:::::=
@PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
+
CCEFEFFFFFFFFFFFE===EEFFDDDFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAAADDFDAAAFFFFFFEEEEEEEEFFFFFFFDDDFFFFFFFFFFFFFFEDDDFFFFFFFFFFEEEFEDDDF===ADFFFEEEFFFFFFDDDDDFFFFEEEF====DEFFFFFFFFFFEBBA@EEBB<<<E@:::AEEEEEEEEEE
@PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
+
IIIIIIIIIIIIFFFFFBBBEEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFBBBBFFFFFFFFFFFFFFEE====BBBEEFFFFFFFFFFFFFFFFFFFFFFFFEEEEEE@@<<===<6@@EEEEEEEEEEE@@@@@@@
@PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
+
26555@@BDDDDDBBBD@@@DDDDDDD@@@BDDDDDDDDDD@@@;;;;DDDDDDDB@@@DDDDDDDDDDDDDD@@@DDDBBBBBDDDDDDDDDDDDDDDDD?;;;?BDDDDDDDDBBBDB<<:::<
'''

expected_fasta_default_options = '''>PC.634_1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
>PC.634_2
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.354_3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.481_4
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
'''

expected_qual_default_options = '''>PC.634_1
40 40 39 39 38 40 40 40 40 40 40 40 37 37 37 37 37 35 35 35 37 37 37 37 37 35 35 35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 28 28 28 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36
36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 35 32 32 32 32 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37 37
37 36 36 31 31 32 32 36 36 36 36 36 36 36 36 36 36 36 28 27 27 27 26 26 26 30 29 30 29 24 24 24 21 15 15 13 13
>PC.634_2
24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37
37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31 31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25
25 25 28
>PC.354_3
34 34 36 37 36 37 37 37 37 37 37 37 37 37 37 37 36 28 28 28 36 36 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 32 32 32 35 35 37 35 32 32 32 37 37 37 37 37 37 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 35 35 35 37 37
37 37 37 37 37 37 37 37 36 36 36 37 36 35 35 35 37 28 28 28 32 35 37 37 37 36 36 36 37 37 37 37 37 37 35 35 35 35 35 37 37 37 37 36 36 36 37 28 28 28 28 35 36 37 37 37 37 37 37 37
37 37 37 36 33 33 32 31 36 36 33 33 27 27 27 36 31 25 25 25 32 36 36 36 36 36 36 36 36 36 36
>PC.481_4
40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
>PC.634_5
17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35
35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33
27 27 25 25 25 27
'''
expected_fasta_full_fasta_headers = '''>PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
'''

expected_qual_full_fasta_headers = '''>PC.634_1 FLP3FBN01ELBSX orig_bc=GCAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=1
40 40 39 39 38 40 40 40 40 40 40 40 37 37 37 37 37 35 35 35 37 37 37 37 37 35 35 35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 28 28 28 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36
36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 35 32 32 32 32 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37 37
37 36 36 31 31 32 32 36 36 36 36 36 36 36 36 36 36 36 28 27 27 27 26 26 26 30 29 30 29 24 24 24 21 15 15 13 13
>PC.634_2 FLP3FBN01EG8AX orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37
37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31 31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25
25 25 28
>PC.354_3 FLP3FBN01EEWKD orig_bc=AGCACGAGCCTA new_bc=AGCACGAGCCTA bc_diffs=0
34 34 36 37 36 37 37 37 37 37 37 37 37 37 37 37 36 28 28 28 36 36 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 32 32 32 35 35 37 35 32 32 32 37 37 37 37 37 37 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 35 35 35 37 37
37 37 37 37 37 37 37 37 36 36 36 37 36 35 35 35 37 28 28 28 32 35 37 37 37 36 36 36 37 37 37 37 37 37 35 35 35 35 35 37 37 37 37 36 36 36 37 28 28 28 28 35 36 37 37 37 37 37 37 37
37 37 37 36 33 33 32 31 36 36 33 33 27 27 27 36 31 25 25 25 32 36 36 36 36 36 36 36 36 36 36
>PC.481_4 FLP3FBN01DEHK3 orig_bc=ACCAGCGACTAG new_bc=ACCAGCGACTAG bc_diffs=0
40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
>PC.634_5 FLP3FBN01DGFYQ orig_bc=ACAGAGTCGGCT new_bc=ACAGAGTCGGCT bc_diffs=0
17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35
35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33
27 27 25 25 25 27
'''

expected_fasta_354_default = '''>PC.354_3
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGATCAGTCTCTTAACTCGGCTATGCATCATTGCCTTGGTAAGCCGTTACCTTACCAACTAGCTAATGCACCGCAGGTCCATCCAAGAGTGATAGCAGAACCATCTTTCAAACTCTAGACATGCGTCTAGTGTTGTTATCCGGTATTAGCATCTGTTTCCAGGTGTTATCCCAGTCTCTTGGG
'''

expected_qual_354_default = '''>PC.354_3
34 34 36 37 36 37 37 37 37 37 37 37 37 37 37 37 36 28 28 28 36 36 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 32 32 32 35 35 37 35 32 32 32 37 37 37 37 37 37 36 36 36 36 36 36 36 36 37 37 37 37 37 37 37 35 35 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 35 35 35 37 37
37 37 37 37 37 37 37 37 36 36 36 37 36 35 35 35 37 28 28 28 32 35 37 37 37 36 36 36 37 37 37 37 37 37 35 35 35 35 35 37 37 37 37 36 36 36 37 28 28 28 28 35 36 37 37 37 37 37 37 37
37 37 37 36 33 33 32 31 36 36 33 33 27 27 27 36 31 25 25 25 32 36 36 36 36 36 36 36 36 36 36
'''

expected_fasta_481_default = '''>PC.481_4
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
'''

expected_qual_481_default = '''>PC.481_4
40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36 36 36 36 36 36 31 31 31 31 31 31 31
'''

expected_fasta_634_default = '''>PC.634_1
CTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCTAGG
>PC.634_2
TTGGACCGTGTCTCAGTTCCAATGTGGGGGCCTTCCTCTCAGAACCCCTATCCATCGAAGGCTTGGTGGGCCGTTACCCCGCCAACAACCTAATGGAACGCATCCCCATCGATGACCGAAGTTCTTTAATAGTTCTACCATGCGGAAGAACTATGCCATCGGGTATTAATCTTTCTTTCGAAAGGCTATCCCCGAGTCATCGGCAGGTTGGATACGTGTTACTCACCCGTGCGCCGGTCGCCA
>PC.634_5
TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATC
'''

expected_qual_634_default = '''>PC.634_1
40 40 39 39 38 40 40 40 40 40 40 40 37 37 37 37 37 35 35 35 37 37 37 37 37 35 35 35 31 31 23 23 23 31 21 21 21 35 35 37 37 37 36 36 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 28 28 28 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 36 36
36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 35 32 32 32 32 35 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 32 32 32 36 37 35 32 32 32 32 32 32 32 32 36 37 37 37
37 36 36 31 31 32 32 36 36 36 36 36 36 36 36 36 36 36 28 27 27 27 26 26 26 30 29 30 29 24 24 24 21 15 15 13 13
>PC.634_2
24 29 33 33 39 39 39 40 39 39 39 40 37 37 37 37 37 37 37 37 37 37 37 32 32 20 20 20 20 20 35 35 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 37 37 37 37 37 36 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 33 28 28 28 28 36 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 36 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 28 28 28 37 28 28 28 37 37
37 37 37 36 36 36 36 36 28 26 26 26 26 28 36 36 36 36 36 36 36 37 38 38 38 38 38 37 37 37 37 37 31 31 31 31 31 31 31 31 31 31 31 31 30 22 22 22 25 25 31 31 31 31 31 31 31 25 25 25
25 25 28
>PC.634_5
17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35 31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35
35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33
27 27 25 25 25 27
'''

if __name__ == '__main__':
    main()
