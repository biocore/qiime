#!/usr/bin/env python

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["William Walters", "Daniel McDonald"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"


from os.path import exists, join
from shutil import rmtree
from re import compile
from tempfile import mkdtemp
from unittest import TestCase, main

import numpy as np
from skbio.util import create_dir
from qiime.extract_barcodes import (extract_barcodes,
                                    process_barcode_single_end_data, process_barcode_paired_end_data,
                                    process_barcode_paired_stitched, process_barcode_in_label,
                                    get_primers)


class FakeOutFile(object):

    def __init__(self, name="test_file"):
        self.data = ""
        self.name = name

    def write(self, s):
        self.data += s


class ExtractBarcodes(TestCase):

    def setUp(self):
        # create the temporary input files that will be used

        self.iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]',
                      'Y': '[CT]', 'S': '[GC]', 'W':
                      '[AT]', 'K': '[GT]', 'M': '[AC]',
                      'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}

        self.output_dir = mkdtemp()
        self.output_dir += '/'

        create_dir(self.output_dir)

    def tearDown(self):
        if exists(self.output_dir):
            rmtree(self.output_dir)

    def test_extract_barcodes_single_end(self):
        """ Extracts barcodes from single-end read """

        fastq_lines =\
            "@HWI-ST830\nAAAATTTTCCCCGGGG\n+\n1234567890ABCDEF".split('\n')

        extract_barcodes(fastq_lines, output_dir=self.output_dir,
                         disable_header_match=True)

        output_bcs_fp = open(join(self.output_dir, "barcodes.fastq"), "U")
        actual_bcs = [line for line in output_bcs_fp]
        expected_bcs = ['@HWI-ST830\n', 'AAAATT\n', '+\n', '123456\n']

        self.assertEqual(actual_bcs, expected_bcs)

        output_reads_fp = open(join(self.output_dir, "reads.fastq"), "U")
        actual_reads = [line for line in output_reads_fp]
        expected_reads = [
            '@HWI-ST830\n',
            'TTCCCCGGGG\n',
            '+\n',
            '7890ABCDEF\n']

        self.assertEqual(actual_reads, expected_reads)

    def test_extract_barcodes_paired_end(self):
        """ Extracts barcodes from paired-end reads """

        fastq1_lines =\
            "@HWI-ST830\nAAAATTTTCCCCGGGG\n+\n1234567890ABCDEF".split('\n')
        fastq2_lines =\
            "@HWI-ST830\nGGGGTTTTAAAACCCC\n+\n1234567890ABCDEF".split('\n')

        extract_barcodes(fastq1=fastq1_lines, fastq2=fastq2_lines,
                         input_type="barcode_paired_end", output_dir=self.output_dir)

        output_bcs_fp = open(join(self.output_dir, "barcodes.fastq"), "U")
        actual_bcs = [line for line in output_bcs_fp]
        expected_bcs =\
            ['@HWI-ST830\n', 'AAAATTGGGGTT\n', '+\n', '123456123456\n']

        self.assertEqual(actual_bcs, expected_bcs)

        # reads 1 output
        output_reads_fp = open(join(self.output_dir, "reads1.fastq"), "U")
        actual_reads = [line for line in output_reads_fp]
        expected_reads = [
            '@HWI-ST830\n',
            'TTCCCCGGGG\n',
            '+\n',
            '7890ABCDEF\n']

        self.assertEqual(actual_reads, expected_reads)
        # reads 2 output
        output_reads_fp = open(join(self.output_dir, "reads2.fastq"), "U")
        actual_reads = [line for line in output_reads_fp]
        expected_reads = [
            '@HWI-ST830\n',
            'TTAAAACCCC\n',
            '+\n',
            '7890ABCDEF\n']

        self.assertEqual(actual_reads, expected_reads)

    def test_extract_barcodes_stitched_reads(self):
        """ Extracts barcodes from ends of a single read """

        fastq_lines =\
            "@HWI-ST830\nAAAATTTTCCCCGGGG\n+\n1234567890ABCDEF\n".split('\n')

        extract_barcodes(fastq_lines, input_type="barcode_paired_stitched",
                         output_dir=self.output_dir, disable_header_match=True)

        output_bcs_fp = open(join(self.output_dir, "barcodes.fastq"), "U")
        actual_bcs = [line for line in output_bcs_fp]
        expected_bcs =\
            ['@HWI-ST830\n', 'AAAATTCCGGGG\n', '+\n', '123456ABCDEF\n']

        self.assertEqual(actual_bcs, expected_bcs)

        output_reads_fp = open(join(self.output_dir, "reads.fastq"), "U")
        actual_reads = [line for line in output_reads_fp]
        expected_reads = ['@HWI-ST830\n', 'TTCC\n', '+\n', '7890\n']

        self.assertEqual(actual_reads, expected_reads)

    def test_extract_barcodes_from_labels(self):
        """ Extracts barcodes from fastq labels """

        fastq_lines =\
            "@HWI-ST830:GTATCT\nAAAATTTTCCCCGGGG\n+\n1234567890ABCDEF".split('\n')

        extract_barcodes(fastq_lines, input_type="barcode_in_label",
                         output_dir=self.output_dir, disable_header_match=True)

        output_bcs_fp = open(join(self.output_dir, "barcodes.fastq"), "U")
        actual_bcs = [line for line in output_bcs_fp]
        expected_bcs =\
            ['@HWI-ST830:GTATCT\n', 'GTATCT\n', '+\n', "''''''\n"]

        self.assertEqual(actual_bcs, expected_bcs)

    def test_process_barcode_single_end_data(self):
        """ Handles fastq lines, parses barcodes """

        fastq_data = ["HWI-ST830", "AAAATTTTCCCCGGGG",
                      np.arange(3, 19, dtype=np.int8)]
        reads_out = FakeOutFile()
        bcs_out = FakeOutFile()

        process_barcode_single_end_data(fastq_data, bcs_out, reads_out,
                                        bc1_len=5, rev_comp_bc1=True)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ["@HWI-ST830", "ATTTT", "+", "('&%$", ""]

        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TTTCCCCGGGG', '+', ')*+,-./0123', '']

        self.assertEqual(actual_reads, expected_reads)

    def test_process_barcode_paired_end_data(self):
        """ Handles paired fastq lines, parses barcodes """

        fastq1_data = ["HWI-ST830", "AAAATTTTCCCCGGGG",
                       np.arange(3, 19, dtype=np.int8)]
        fastq2_data = ["HWI-ST830", "TCCCCGGGG", np.arange(3, 12, dtype=np.int8)]
        reads1_out = FakeOutFile()
        reads2_out = FakeOutFile()
        bcs_out = FakeOutFile()

        process_barcode_paired_end_data(fastq1_data, fastq2_data,
                                        bcs_out, reads1_out, reads2_out, bc1_len=5, bc2_len=3,
                                        rev_comp_bc1=True, rev_comp_bc2=True)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'ATTTTGGA', '+', "('&%$&%$", '']

        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TTTCCCCGGGG', '+', ')*+,-./0123', '']

        self.assertEqual(actual_reads, expected_reads)

        actual_reads = reads2_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'CCGGGG', '+', "'()*+,", '']

        self.assertEqual(actual_reads, expected_reads)

    def test_process_barcode_paired_end_data_orientation_no_match(self):
        """ Handles paired fastq lines, parses barcodes, orients reads """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        fastq2_data = ["HWI-ST830", "GGTTCCAA", np.arange(3, 11, dtype=np.int8)]
        reads1_out = FakeOutFile()
        reads2_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'AYA']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'ATA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()
        fastq2_out_not_oriented = FakeOutFile()

        # With no matches, should write to the not_oriented files, and keep
        # in the same order of file 1 and file 2
        process_barcode_paired_end_data(fastq1_data, fastq2_data,
                                        bcs_out, reads1_out, reads2_out, bc1_len=5, bc2_len=3,
                                        rev_comp_bc1=False, rev_comp_bc2=False,
                                        attempt_read_orientation=True, forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq1_out_not_oriented=fastq1_out_not_oriented,
                                        fastq2_out_not_oriented=fastq2_out_not_oriented)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads, expected_reads)

        actual_reads = reads2_out.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'ATCGAGGT', '+', "$%&'($%&", '']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCGATCGATCGATCG', '+',
                          ')*+,-./01234567', '']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

        actual_reads_not_oriented = fastq2_out_not_oriented.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCCAA', '+', "'()*+", '']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_end_data_orientation_forward_match(self):
        """ Handles paired fastq lines, parses barcodes, orients reads """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        fastq2_data = ["HWI-ST830", "GGTTCCAA", np.arange(3, 11, dtype=np.int8)]
        reads1_out = FakeOutFile()
        reads2_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'GATCGA']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'ATA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()
        fastq2_out_not_oriented = FakeOutFile()

        # With a match to the forward primer, should parse out primers in
        # the given order of read 1 and read 2.
        process_barcode_paired_end_data(fastq1_data, fastq2_data,
                                        bcs_out, reads1_out, reads2_out, bc1_len=5, bc2_len=3,
                                        rev_comp_bc1=False, rev_comp_bc2=False,
                                        attempt_read_orientation=True, forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq1_out_not_oriented=fastq1_out_not_oriented,
                                        fastq2_out_not_oriented=fastq2_out_not_oriented)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'ATCGAGGT', '+', "$%&'($%&", '']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCGATCGATCGATCG', '+',
                          ')*+,-./01234567', '']
        self.assertEqual(actual_reads, expected_reads)

        actual_reads = reads2_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCCAA', '+', "'()*+", '']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

        actual_reads_not_oriented = fastq2_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_end_data_orientation_forward_in_read2(
            self):
        """ Handles paired fastq lines, parses barcodes, orients reads """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        fastq2_data = ["HWI-ST830", "GGTTCCAA", np.arange(3, 11, dtype=np.int8)]
        reads1_out = FakeOutFile()
        reads2_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'TTCCA']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'ATA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()
        fastq2_out_not_oriented = FakeOutFile()

        # With a forward primer match in read 2, should reverse read order
        process_barcode_paired_end_data(fastq1_data, fastq2_data,
                                        bcs_out, reads1_out, reads2_out, bc1_len=5, bc2_len=3,
                                        rev_comp_bc1=False, rev_comp_bc2=False,
                                        attempt_read_orientation=True, forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq1_out_not_oriented=fastq1_out_not_oriented,
                                        fastq2_out_not_oriented=fastq2_out_not_oriented)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'GGTTCATC', '+', "$%&'($%&", '']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'CAA', '+', ')*+', '']
        self.assertEqual(actual_reads, expected_reads)

        actual_reads = reads2_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'GATCGATCGATCGATCG', '+',
                          "'()*+,-./01234567", '']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

        actual_reads_not_oriented = fastq2_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_end_data_orientation_reverse_in_read1(
            self):
        """ Handles paired fastq lines, parses barcodes, orients reads """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        fastq2_data = ["HWI-ST830", "GGTTCCAA", np.arange(3, 11, dtype=np.int8)]
        reads1_out = FakeOutFile()
        reads2_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'TTTTT']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'CGATCGA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()
        fastq2_out_not_oriented = FakeOutFile()

        # With a forward primer match in read 2, should reverse read order
        process_barcode_paired_end_data(fastq1_data, fastq2_data,
                                        bcs_out, reads1_out, reads2_out, bc1_len=5, bc2_len=3,
                                        rev_comp_bc1=False, rev_comp_bc2=False,
                                        attempt_read_orientation=True, forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq1_out_not_oriented=fastq1_out_not_oriented,
                                        fastq2_out_not_oriented=fastq2_out_not_oriented)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'GGTTCATC', '+', "$%&'($%&", '']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'CAA', '+', ')*+', '']
        self.assertEqual(actual_reads, expected_reads)

        actual_reads = reads2_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'GATCGATCGATCGATCG', '+',
                          "'()*+,-./01234567", '']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

        actual_reads_not_oriented = fastq2_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_end_data_orientation_rev_in_read2(self):
        """ Handles paired fastq lines, parses barcodes, orients reads """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        fastq2_data = ["HWI-ST830", "GGTTCCAA", np.arange(3, 11, dtype=np.int8)]
        reads1_out = FakeOutFile()
        reads2_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'TTTTTT']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'TCCAA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()
        fastq2_out_not_oriented = FakeOutFile()

        # With a reverse primer in read 2, should write in current order.
        process_barcode_paired_end_data(fastq1_data, fastq2_data,
                                        bcs_out, reads1_out, reads2_out, bc1_len=5, bc2_len=3,
                                        rev_comp_bc1=False, rev_comp_bc2=False,
                                        attempt_read_orientation=True, forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq1_out_not_oriented=fastq1_out_not_oriented,
                                        fastq2_out_not_oriented=fastq2_out_not_oriented)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'ATCGAGGT', '+', "$%&'($%&", '']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCGATCGATCGATCG', '+',
                          ')*+,-./01234567', '']
        self.assertEqual(actual_reads, expected_reads)

        actual_reads = reads2_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCCAA', '+', "'()*+", '']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

        actual_reads_not_oriented = fastq2_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_stitched(self):
        """ Handles stitched barcode data, parses barcodes from ends """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        reads1_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'ATA']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'ATA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()

        # With no matches, should write to the not_oriented files, and keep
        # in the same order of output file
        process_barcode_paired_stitched(fastq1_data,
                                        bcs_out, reads1_out, bc1_len=3, bc2_len=4,
                                        rev_comp_bc1=True, rev_comp_bc2=True,
                                        attempt_read_orientation=True,
                                        forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq_out_not_oriented=fastq1_out_not_oriented,
                                        switch_bc_order=True)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'CGATGAT', '+', '7654&%$', '']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads =\
            ['@HWI-ST830', 'GATCGATCGATCG', '+', "'()*+,-./0123", '']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_stitched_forward_primer_match(self):
        """ Handles stitched barcode data, parses barcodes from ends """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        reads1_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'GATCGA']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'ATA']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()

        # With forward primer match, should write in order of read 1, read 2
        process_barcode_paired_stitched(fastq1_data,
                                        bcs_out, reads1_out, bc1_len=3, bc2_len=4,
                                        rev_comp_bc1=True, rev_comp_bc2=True,
                                        attempt_read_orientation=True,
                                        forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq_out_not_oriented=fastq1_out_not_oriented,
                                        switch_bc_order=True)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'CGATGAT', '+', '7654&%$', '']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'GATCGATCGATCG', '+',
                          "'()*+,-./0123", '']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_paired_stitched_reverse_primer_match(self):
        """ Handles stitched barcode data, parses barcodes from ends """

        fastq1_data = ["HWI-ST830", "ATCGATCGATCGATCGATCG",
                       np.arange(3, 23, dtype=np.int8)]
        reads1_out = FakeOutFile()
        bcs_out = FakeOutFile()
        forward_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'AAAAAA']))]
        reverse_primers = [compile(''.join([self.iupac[symbol] for
                                            symbol in 'GATCG']))]
        output_bc_not_oriented = FakeOutFile()
        fastq1_out_not_oriented = FakeOutFile()

        # With reverse primer match, should write in order of read2, read 1
        process_barcode_paired_stitched(fastq1_data,
                                        bcs_out, reads1_out, bc1_len=3, bc2_len=4,
                                        rev_comp_bc1=True, rev_comp_bc2=False,
                                        attempt_read_orientation=True,
                                        forward_primers=forward_primers,
                                        reverse_primers=reverse_primers,
                                        output_bc_not_oriented=output_bc_not_oriented,
                                        fastq_out_not_oriented=fastq1_out_not_oriented,
                                        switch_bc_order=False)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830', 'TCGCGAT', '+', "567'&%$", '']
        self.assertEqual(actual_bcs, expected_bcs)

        actual_reads = reads1_out.data.split('\n')
        expected_reads = ['@HWI-ST830', 'TCGATCGATCGAT', '+',
                          '43210/.-,+*)(', '']
        self.assertEqual(actual_reads, expected_reads)

        actual_bcs_not_oriented = output_bc_not_oriented.data.split('\n')
        expected_bcs = ['']
        self.assertEqual(actual_bcs_not_oriented, expected_bcs)

        actual_reads_not_oriented = fastq1_out_not_oriented.data.split('\n')
        expected_reads = ['']
        self.assertEqual(actual_reads_not_oriented, expected_reads)

    def test_process_barcode_in_label(self):
        """ Handles label barcodes from one or two fastq labels """

        fastq1_data = [
            "HWI-ST830:ATCG",
            "AAAATTTTCCCCGGGG",
            np.arange(3, 19, dtype=np.int8)]
        fastq2_data = ["HWI-ST830:GGGG", "TCCCCGGGG",
                       np.arange(3, 12, dtype=np.int8)]
        bcs_out = FakeOutFile()

        process_barcode_in_label(fastq1_data, fastq2_data, bcs_out,
                                 bc1_len=4, bc2_len=3, rev_comp_bc1=True, rev_comp_bc2=True)

        actual_bcs = bcs_out.data.split('\n')
        expected_bcs = ['@HWI-ST830:ATCG', 'CGATCCC', '+', "'''''''", '']
        self.assertEqual(actual_bcs, expected_bcs)

    def test_get_primers(self):
        """ Get primer regular expression generators out of mapping data """

        # Raise error if ReversePrimer not supplied
        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'Description']
        mapping_data = [['s1', 'ATCG', 'TTGGCC,TTGGWC', 'ATRCCTA']]
        self.assertRaises(IndexError, get_primers, header, mapping_data)

        header = ['SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                  'ReversePrimer', 'Description']
        forward_primers, reverse_primers = get_primers(header, mapping_data)

        forward_primers = set([seq.pattern for seq in forward_primers])
        reverse_primers = set([seq.pattern for seq in reverse_primers])

        expected_forward_primers = set(['TTGGCC', 'TAGG[CT]AT', 'TTGG[AT]C'])
        expected_reverse_primers = set(['GGCCAA', 'AT[AG]CCTA', 'G[AT]CCAA'])

        self.assertEqual(forward_primers, expected_forward_primers)
        self.assertEqual(reverse_primers, expected_reverse_primers)


if __name__ == '__main__':
    main()
