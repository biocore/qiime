#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the remote.py module."""

from cogent.util.unit_test import TestCase, main
from qiime.remote import (_extract_spreadsheet_key_from_url,
                          load_google_spreadsheet_mapping_file,
                          RemoteMappingFileError)

class RemoteTests(TestCase):
    """Tests for the remote.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.url = url
        self.spreadsheet_key = '0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc'
        self.worksheet_name = 'Fasting_Map'
        self.exp_mapping_lines = exp_mapping_lines

    def test_load_google_spreadsheet_mapping_file(self):
        """Test correctly retrieves a remote mapping file."""
        # If we have an active Internet connection, try the test. If not, at
        # least make sure an appropriate error is thrown.
        try:
            obs = load_google_spreadsheet_mapping_file(self.spreadsheet_key,
                                                       worksheet_name=None)
        except RemoteMappingFileError:
            pass
        else:
            self.assertEqual(obs, self.exp_mapping_lines)

    def test_extract_spreadsheet_key_from_url(self):
        """Test correctly extracts a key from a URL."""
        # Pass a URL.
        obs = _extract_spreadsheet_key_from_url(self.url)
        self.assertEqual(obs, self.spreadsheet_key)

        # Pass a key.
        obs = _extract_spreadsheet_key_from_url(self.spreadsheet_key)
        self.assertEqual(obs, self.spreadsheet_key)


url = 'https://docs.google.com/spreadsheet/ccc?key=0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc#gid=1'

exp_mapping_lines = """#SampleID	BarcodeSequence	LinkerPrimerSequence	Treatment	DOB	Description
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._354
PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	Control	20061218	Control_mouse_I.D._355
PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	Control	20061126	Control_mouse_I.D._356
PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	Control	20070314	Control_mouse_I.D._481
PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	Control	20071210	Control_mouse_I.D._593
PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	Fast	20071112	Fasting_mouse_I.D._607
PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._634
PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._635
PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	Fast	20080116	Fasting_mouse_I.D._636
"""


if __name__ == "__main__":
    main()
