#!/usr/bin/env python

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import isdir, isfile, basename
from shutil import rmtree
from os import close
from tempfile import mkstemp

from numpy import array
from unittest import TestCase, main
from skbio.util import remove_files

from qiime.util import create_dir
from qiime.truncate_fasta_qual_files import parse_fasta_file,\
    verify_equivalency, truncate_fasta_qual, truncate_seqs, get_output_filepaths,\
    write_trunc_fasta, write_trunc_qual


class TruncateFastaQualFilesTests(TestCase):

    """ Unit tests for the truncate_fasta_qual_files.py module """

    def setUp(self):
        # create the temporary input files that will be used

        self._files_to_remove = []

        # create the temporary input files that will be used with the
        # filter_fasta function
        fd, self.fasta_fp = mkstemp(prefix='fasta_seqs_',
                                    suffix='.fasta')
        close(fd)
        seq_file = open(self.fasta_fp, 'w')
        seq_file.write(fasta_seqs)
        seq_file.close()

        fd, self.qual_fp = mkstemp(prefix='qual_scores_',
                                   suffix='.qual')
        close(fd)
        seq_file = open(self.qual_fp, 'w')
        seq_file.write(qual_scores)
        seq_file.close()

        fd, self.diff_name_fasta_fp = mkstemp(prefix='fasta_seqs_diff_name_',
                                              suffix='.fasta')
        close(fd)
        seq_file = open(self.diff_name_fasta_fp, 'w')
        seq_file.write(diff_name_fasta_seqs)
        seq_file.close()

        fd, self.diff_len_fasta_fp = mkstemp(prefix='fasta_seqs_diff_len_',
                                             suffix='.fasta')
        close(fd)
        seq_file = open(self.diff_len_fasta_fp, 'w')
        seq_file.write(diff_len_fasta_seqs)
        seq_file.close()

        self._files_to_remove =\
            [self.fasta_fp, self.qual_fp, self.diff_name_fasta_fp,
             self.diff_len_fasta_fp]

    def tearDown(self):
        if self._files_to_remove:
            remove_files(self._files_to_remove)
        if isdir('/tmp/truncate_fasta_qual_test/'):
            rmtree('/tmp/truncate_fasta_qual_test/')

    def test_parse_fasta_file(self):
        """ Properly returns dict from fasta data """

        fasta_data = ['>seq1 SAMPLE1', 'AAACGT', '>seq2', 'ACGGT']

        expected_fasta = {'seq1': 'AAACGT', 'seq2': 'ACGGT'}

        expected_order = ['seq1 SAMPLE1', 'seq2']

        actual_fasta, actual_order = parse_fasta_file(fasta_data)

        self.assertEqual(actual_fasta, expected_fasta)

        self.assertEqual(actual_order, expected_order)

    def test_verify_equivalency(self):
        """ verifies identical labels, base counts between qual and fasta """

        fasta_data = {'seq1': 'AATC', 'seq2': 'GGAT'}

        qual_data = {'seq1': ["40", "36", "35", "18"],
                     'seq2': ["39", "18", "22", "23"]}

        # Should not raise any errors
        verify_equivalency(fasta_data, qual_data)

        # Make number of sequences unequal, and should raise error.
        fasta_data = {'seq1': 'AATC', 'seq2': 'GGAT', 'seq3': 'ACTG'}

        qual_data = {'seq1': ["40", "36", "35", "18"],
                     'seq2': ["39", "18", "22", "23"]}

        self.assertRaises(
            ValueError,
            verify_equivalency,
            fasta_data,
            qual_data)

        # Should raise error if lens of values different
        fasta_data = {'seq1': 'AATC', 'seq2': 'GGAT'}

        qual_data = {'seq1': ["50", "40", "36", "35", "18"],
                     'seq2': ["39", "18", "22", "23"]}

        self.assertRaises(
            ValueError,
            verify_equivalency,
            fasta_data,
            qual_data)

    def test_truncate_fasta_qual(self):
        """ Test for overall module functionality """

        base_pos = 80
        output_dir = '/tmp/truncate_fasta_qual_test/'

        create_dir(output_dir)

        truncate_fasta_qual(self.fasta_fp, self.qual_fp, output_dir, base_pos)

        actual_trunc_fasta_fp = output_dir +\
            basename(self.fasta_fp).replace(".fasta", "_filtered.fasta")

        actual_trunc_fasta_fp = open(actual_trunc_fasta_fp, "U")

        actual_trunc_fasta = [line.strip() for line in actual_trunc_fasta_fp]

        self.assertEqual(actual_trunc_fasta, expected_fasta_seqs)

        actual_trunc_qual_fp = output_dir +\
            basename(self.qual_fp).replace(".qual", "_filtered.qual")

        actual_trunc_qual_fp = open(actual_trunc_qual_fp, "U")

        actual_trunc_qual = [line.strip() for line in actual_trunc_qual_fp]

        self.assertEqual(actual_trunc_qual, expected_qual_scores)

    def test_truncate_seqs(self):
        """ Truncates seqs, scores correctly """

        base_pos = 5

        fasta_seqs = {'seq1': 'GAAATCAAGAATAC',
                      'seq2': 'ATAAACAAGAT'}
        qual_scores = {'seq1': array(map(str, [20, 10, 15, 25, 24, 25, 27])),
                       'seq2': array(map(str, [22, 21, 15, 12, 22, 25, 27, 28]))}

        expected_fasta = {'seq1': 'GAAAT',
                          'seq2': 'ATAAA'}

        expected_qual = {'seq1': map(str, array([20, 10, 15, 25, 24])),
                         'seq2': map(str, array([22, 21, 15, 12, 22]))}

        actual_fasta_seqs, actual_qual_scores =\
            truncate_seqs(fasta_seqs, qual_scores, base_pos)

        self.assertDictEqual(actual_fasta_seqs, expected_fasta)
        self.assertItemsEqual(expected_qual, actual_qual_scores)

    def test_get_output_filepaths(self):
        """ Generates output filepaths for fasta, qual files correctly """

        output_dir = "."

        fasta_fp = "seqs.fna"

        qual_fp = "seqs.qual"

        expected_fasta_fp = "./seqs_filtered.fasta"
        expected_qual_fp = "./seqs_filtered.qual"

        actual_fasta_fp, actual_qual_fp =\
            get_output_filepaths(output_dir, fasta_fp, qual_fp)

        self.assertEqual(actual_fasta_fp, expected_fasta_fp)
        self.assertEqual(actual_qual_fp, expected_qual_fp)

        # Test for relative paths
        output_dir = "test/"

        fasta_fp = "../seqs.fna"

        qual_fp = "quality_scores/seqs.qual"

        expected_fasta_fp = "test/seqs_filtered.fasta"
        expected_qual_fp = "test/seqs_filtered.qual"

        actual_fasta_fp, actual_qual_fp =\
            get_output_filepaths(output_dir, fasta_fp, qual_fp)

        self.assertEqual(actual_fasta_fp, expected_fasta_fp)
        self.assertEqual(actual_qual_fp, expected_qual_fp)

    def test_write_trunc_fasta(self):
        """ writes fasta seq to output filepath correctly """

        seq_order = ['seq1', 'seq2', 'seq3']

        seqs = {'seq1': 'ATCG', 'seq3': 'ACCC', 'seq2': 'GGACC'}

        output_dir = '/tmp/truncate_fasta_qual_test/'

        create_dir(output_dir)

        fasta_out_fp = output_dir + 'seqs_filtered.fna'

        write_trunc_fasta(seqs, fasta_out_fp, seq_order)

        expected_seqs = ['>seq1', 'ATCG', '>seq2', 'GGACC', '>seq3', 'ACCC']

        actual_fasta = open(fasta_out_fp, "U")

        actual_fasta = [line.strip() for line in actual_fasta]

        self.assertEqual(actual_fasta, expected_seqs)

    def test_write_trunc_qual(self):
        """ writes truncated qual scores out in correct format """

        seq_order = ['seq1', 'seq2', 'seq3']

        output_dir = '/tmp/truncate_fasta_qual_test/'

        create_dir(output_dir)

        qual_out_fp = output_dir + 'seqs_filtered.qual'

        write_trunc_qual(trunc_qual_scores, qual_out_fp, seq_order)

        # Needs to correctly insert newline after every 60 base calls
        expected_scores = expected_qual_scores

        actual_qual = open(qual_out_fp, "U")

        actual_qual = [line.strip() for line in actual_qual]

        self.assertEqual(actual_qual, expected_scores)


# Long strings at the end for better readability
fasta_seqs = """>seq1
ACCAGCGACTAGCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>seq2
ACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>seq3
ACGGTGAGTGTCCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG"""

expected_fasta_seqs = """>seq1
ACCAGCGACTAGCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCG
>seq2
ACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCG
>seq3
ACGGTGAGTGTCCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACC""".split('\n')

qual_scores = """>seq1
36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 40 40 40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 33 33 33 33 37 37 37 37 37 37 37 37 37 37 37 37
37 37 36 36 28 28 28 28 33 33 33 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 36 36 36 36 36 36 31 31 27 27 28 28 28 27 21 31 31 36 36 36 36 36 36
36 36 36 36 36 31 31 31 31 31 31 31
>seq2
35 35 35 35 33 31 31 31 33 35 35 35 35 35 35 35 35 35 35 35 35 35 23 20 20 31 31 33 33 33 35 23 17 17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35
31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26 26 35 35 35 35 35 35 35 33 31 31 31 35 35 35 35 35 35 35 35 35 35 35 35 35 35 31 31 31 35 35 35 33 33 33 33 33 35 35 35
35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 26 26 26 30 33 35 35 35 35 35 35 35 35 33 33 33 35 33 27 27 25 25 25 27 14 14 14 14 14 25 25 34 34 35 35 35 32 33 33 32 35 35 32 25 25
15 20 20 20 28 35 33 33 33 33 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 29 24 24 24 29 35 35 35 35 33 33 31 31 34 34 34 34 34 34 31 20 20 20 20 20 31 34 31 31 31 31 32 31 31 33
34 25 25 20 20 18 25 28 28 22 20 22 28 28 28 30 30 29 29 29 30 25 25 25 29 29 26 26 25
>seq3
32 32 32 32 35 35 35 35 35 35 35 35 35 35 35 35 35 35 38 38 39 39 32 32 32 35 35 35 35 35 34 31 21 21 25 35 32 25 25 25 32 35 35 37 39 35 35 35 35 35 35 35 35 35 35 35 35 35 32 32
32 32 32 35 32 32 32 32 35 35 35 35 35 35 35 35 35 35 32 32 32 32 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 34 34 26 26 26 26 32 35 35 35 35 35
35 34 34 34 34 34 34 34 35 35 35 35 35 35 35 35 35 26 26 26 26 35 35 35 35 35 35 35 35 34 34 34 34 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 30 22 24 21 21 21 30 35 35 35 35
35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 35 29 29 26 34 35 27 27 27 27 35 35 35 35 35 35 35 32 32 32 32 32 35 35 35 35 35 35 35 35
35 35 31 32 32 25 28 25 25 25 25 30 30 30 30 30 30 30 30 28 22 22 22 28 28 30 25 22 22 22 30 30"""

expected_qual_scores = """>seq1
36 36 36 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 40 40 40 40 40 40 40 40 40 40 38 38 39 40 40 40 40 40 40 40 40 40 40 40 40 40 40 37 37 37 37 37 33 33 33 36 36 37 37 37 37 37
37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37 37
>seq2
35 35 35 35 33 31 31 31 33 35 35 35 35 35 35 35 35 35 35 35 35 35 23 20 20 31 31 33 33 33 35 23 17 17 21 20 20 20 31 31 33 35 35 35 35 35 33 33 33 35 31 31 31 35 35 35 35 35 35 35
31 31 31 33 35 35 35 35 35 35 35 35 35 35 31 31 31 26 26 26
>seq3
32 32 32 32 35 35 35 35 35 35 35 35 35 35 35 35 35 35 38 38 39 39 32 32 32 35 35 35 35 35 34 31 21 21 25 35 32 25 25 25 32 35 35 37 39 35 35 35 35 35 35 35 35 35 35 35 35 35 32 32
32 32 32 35 32 32 32 32 35 35 35 35 35 35 35 35 35 35 32 32""".split('\n')

diff_name_fasta_seqs = """>DifferentName
ACCAGCGACTAGCATGCTGCCTCCCGTAGGAGTCTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCAACCTCTCAGTCCGGCTACTGATCGTCGACTTGGTGAGCCGTTACCTCACCAACTATCTAATCAGACGCGAGCCCATCTTTCAGCGGATTGCTCCTTTGGTATTCCGGCGATGCCGCCAAAATCATTATGCGGTATTAGCAGTCGTTTCCAACTGTTGTCCCCCTCTGAAAGGCAGGTTGCTCACG
>seq2
ACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>seq3
ACGGTGAGTGTCCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG"""


diff_len_fasta_seqs = """>seq1
ACCAGCGACTAGCATGCTGCCT
>seq2
ACAGAGTCGGCTCATGCTGCCTCCCGTAGGAGTTTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTCACCCTCTCAGGCCGGCTACTGATCGTCGCTTTGGTAGGCCGTTACCCTGCCAACTGGCTAATCAGACGCGGGTCCATCTCACACCGATTAATCTTTTTCCAACCAGAGCATGCGCCCCTGTTGGCTTATGCGGTATTAGCGGTCGTTTCCAACTGTTATCCCCCTGTGTGAGGCAGGTTACCCACGCGTTACTCACCCGTCCG
>seq3
ACGGTGAGTGTCCATGCTGCCTCCCGTAGGAGTTTGGACCGTGTCTCAGTTCCAATGTGGGGGACCTTCCTCTCAGAACCCCTATCCATCGAAGACTAGGTGGGCCGTTACCCCGCCTACTATCTAATGGAACGCATCCCCATCTTATACCGGTAAACCTTTAATCATGAGAAAATGCTCACTCATGATACCATCTTGTATTAATCTCCCTTTCAGAAGGCTATCCAAGAGTATAAGGCAGGTTGGATACGCGTTACTCACCCGTGCGCCGG"""

trunc_qual_scores = {
    'seq1': array(map(str, [36, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 38, 38, 39, 40, 40, 40, 40, 40, 40,
                  40, 40, 40, 40, 40, 40, 40, 40, 37, 37, 37, 37, 37, 33, 33, 33, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37])),
    'seq2': array(map(str, [35, 35, 35, 35, 33, 31, 31, 31, 33, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 23, 20, 20, 31, 31, 33, 33, 33, 35, 23, 17, 17, 21, 20, 20,
                            20, 31, 31, 33, 35, 35, 35, 35, 35, 33, 33, 33, 35, 31, 31, 31, 35, 35, 35, 35, 35, 35, 35, 31, 31, 31, 33, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 31, 31, 31, 26, 26, 26])),
    'seq3': array(map(str, [32, 32, 32, 32, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 38, 38, 39, 39, 32, 32, 32, 35, 35, 35, 35, 35, 34, 31, 21, 21, 25, 35, 32, 25, 25, 25, 32, 35, 35, 37, 39, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 32, 32, 32, 32, 32, 35, 32, 32, 32, 32, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 32, 32]))}

if __name__ == "__main__":
    main()
