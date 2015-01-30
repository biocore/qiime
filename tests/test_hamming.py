#!/usr/bin/env python
# test_hamming.py
"""
Tests for standalone Hamming decoding function

Author: Micah Hamady (hamady@colorado.edu)

"""
__author__ = "Micah Hamady"
__copyright__ = "Copyright 2011, The QIIME Project"  # consider project name
__credits__ = ["Micah Hamady", "Rob Knight"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from unittest import TestCase, main
from qiime.hamming import decode_barcode_8


class GeneralSetUp(TestCase):

    def setUp(self):
        """General Setup"""

        # valid code words
        self.valid_bc_1 = "AACCATGC"
        self.valid_bc_2 = "TCGTAGCA"
        self.valid_bc_3 = "ACACCTCT"
        self.valid_bc_4 = "CTTCCTAG"
        self.valid_bc_5 = "GGTAGCTT"

        # single error reference
        self.single_error_ref = "AACCATGC"

        # A->C
        self.single_error_1 = "ACCCATGC"

        # A->G
        self.single_error_2 = "AACCGTGC"

        # C->A
        self.single_error_3 = "AAACATGC"

        # C->T
        self.single_error_4 = "AACCATGT"

        # T->C
        self.single_error_5 = "AACCACGC"

        # T->G
        self.single_error_6 = "AACCAGGC"

        # G->T
        self.single_error_7 = "AACCATTC"

        # G->A
        self.single_error_8 = "AACCATAC"

        # double error reference
        self.double_error_ref = "AACCATGC"

        # A->T
        self.double_error_1 = "ATCCATGC"

        # T->A
        self.double_error_2 = "AACCAAGC"

        # C->G
        self.double_error_3 = "AACGATGC"

        # aG->C
        self.double_error_4 = "AACCATCC"


class StandaloneHammingTests(GeneralSetUp):

    """Tests for the standalone_hamming functions"""

    def test_decode_barcode_8_ok(self):
        """ Should decode valid codewords w/o error  """
        self.assertEqual(decode_barcode_8(self.valid_bc_1),
                         (self.valid_bc_1, 0))
        self.assertEqual(decode_barcode_8(self.valid_bc_2),
                         (self.valid_bc_2, 0))
        self.assertEqual(decode_barcode_8(self.valid_bc_3),
                         (self.valid_bc_3, 0))
        self.assertEqual(decode_barcode_8(self.valid_bc_4),
                         (self.valid_bc_4, 0))
        self.assertEqual(decode_barcode_8(self.valid_bc_5),
                         (self.valid_bc_5, 0))

    def test_decode_barcode_8_one_error(self):
        """ Should correct single bit errors w/o error """
        self.assertEqual(decode_barcode_8(self.single_error_1),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_2),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_3),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_4),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_5),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_6),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_7),
                         (self.single_error_ref, 0.5))
        self.assertEqual(decode_barcode_8(self.single_error_8),
                         (self.single_error_ref, 0.5))

    def test_decode_barcode_8_two_error(self):
        """ Should raise error when double error detected """
        self.assertEqual(decode_barcode_8(self.double_error_1), (None, 1))
        self.assertEqual(decode_barcode_8(self.double_error_1), (None, 1))
        self.assertEqual(decode_barcode_8(self.double_error_1), (None, 1))
        self.assertEqual(decode_barcode_8(self.double_error_1), (None, 1))

if __name__ == '__main__':
    main()
