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

from socket import gaierror
from cogent.util.unit_test import TestCase, main
from gdata.spreadsheet.service import SpreadsheetsService
from qiime.remote import (_get_cleaned_headers,
                          _get_spreadsheet_headers,
                          _export_mapping_file,
                          _extract_spreadsheet_key_from_url,
                          load_google_spreadsheet_mapping_file,
                          RemoteMappingFileConnectionError,
                          RemoteMappingFileError)

class RemoteTests(TestCase):
    """Tests for the remote.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.url1 = url1
        self.url2 = url2
        self.url3 = url3
        self.spreadsheet_key = '0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc'
        self.worksheet_name = 'Fasting_Map'
        self.worksheet_id = 'od7'
        self.exp_mapping_lines = exp_mapping_lines
        
        # Try to connect if we have an active Internet connection. If not, we
        # can't run most of these tests...
        self.client = SpreadsheetsService()
        try:
            self.client.GetWorksheetsFeed(self.spreadsheet_key,
                                          visibility='public',
                                          projection='basic')
        except gaierror:
            self.client = None

    def test_load_google_spreadsheet_mapping_file(self):
        """Test correctly retrieves a remote mapping file."""
        # If we have an active Internet connection, try the test. If not, at
        # least make sure an appropriate error is thrown.
        try:
            obs = load_google_spreadsheet_mapping_file(self.spreadsheet_key,
                                                       worksheet_name=None)
        except RemoteMappingFileConnectionError:
            pass
        else:
            self.assertEqual(obs, self.exp_mapping_lines)

        # Test naming a worksheet.
        try:
            obs = load_google_spreadsheet_mapping_file(self.spreadsheet_key,
                    worksheet_name='Fasting_Map')
        except RemoteMappingFileConnectionError:
            pass
        else:
            self.assertEqual(obs, self.exp_mapping_lines)

    def test_load_google_spreadsheet_mapping_file_invalid_input(self):
        """Test correctly raises errors on various bad inputs."""
        # Bad worksheet name.
        try:
            self.assertRaises(RemoteMappingFileError,
                    load_google_spreadsheet_mapping_file, self.spreadsheet_key,
                    worksheet_name='foo')
        except RemoteMappingFileConnectionError:
            pass

    def test_extract_spreadsheet_key_from_url(self):
        """Test correctly extracts a key from a URL."""
        # Pass various URLs with different key/value combos.
        obs = _extract_spreadsheet_key_from_url(self.url1)
        self.assertEqual(obs, self.spreadsheet_key)

        obs = _extract_spreadsheet_key_from_url(self.url2)
        self.assertEqual(obs, self.spreadsheet_key)

        obs = _extract_spreadsheet_key_from_url(self.url3)
        self.assertEqual(obs, self.spreadsheet_key)

        # Pass a key directly.
        obs = _extract_spreadsheet_key_from_url(self.spreadsheet_key)
        self.assertEqual(obs, self.spreadsheet_key)

        # Pass 'key=<key>'.
        obs = _extract_spreadsheet_key_from_url('key=' + self.spreadsheet_key)
        self.assertEqual(obs, self.spreadsheet_key)

    def test_get_spreadsheet_headers(self):
        """Test reading the header line from a spreadsheet."""
        if self.client:
            exp = ['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                   'Treatment', 'DOB', 'Description']
            obs = _get_spreadsheet_headers(self.client, self.spreadsheet_key,
                                           self.worksheet_id)
            self.assertEqual(obs, exp)

    def test_export_mapping_file(self):
        """Test exporting mapping file from spreadsheet."""
        if self.client:
            exp = [['#SampleID', 'DOB'],
                ['#Example mapping file for the QIIME analysis package.  '
                 'These 9 samples are from a study of the effects of exercise '
                 'and diet on mouse cardiac physiology (Crawford, et al, '
                 'PNAS, 2009).'], ['PC.354', '20061218'],
                ['PC.355', '20061218'], ['PC.356', '20061126'],
                ['PC.481', '20070314'], ['PC.593', '20071210'],
                ['PC.607', '20071112'], ['PC.634', '20080116'],
                ['PC.635', '20080116'], ['PC.636', '20080116']]
            obs = _export_mapping_file(self.client, self.spreadsheet_key,
                                       self.worksheet_id, ['#SampleID', 'DOB'])
            self.assertEqual(obs, exp)

    def test_export_mapping_file_invalid_input(self):
        """Test exporting mapping file with bad input raises errors."""
        if self.client:
            # Nonexisting header.
            self.assertRaises(RemoteMappingFileError, _export_mapping_file,
                    self.client, self.spreadsheet_key, self.worksheet_id,
                    ['#SampleID', 'Foo'])

    def test_get_cleaned_headers(self):
        """Test correctly converts headers to Google's representation."""
        # Some duplicates.
        exp = ['foo', 'foo_2', 'foo_3', 'foo_4', 'fooo', 'foo_5', 'foo_6']
        obs = _get_cleaned_headers(
                ['foo', 'Foo', 'FOO', 'F_oO', 'F:Oo_o', '#Foo', 'f O O#'])
        self.assertEqual(obs, exp)

        # All unique.
        exp = ['foo', 'bar']
        obs = _get_cleaned_headers(['Fo#o', 'bar'])
        self.assertEqual(obs, exp)

        # Header consisting of only special characters and header that is
        # blank.
        self.assertRaises(RemoteMappingFileError, _get_cleaned_headers,
                          ['Foo', '___', 'BAR'])
        self.assertRaises(RemoteMappingFileError, _get_cleaned_headers,
                          ['Foo', '', 'BAR'])


url1 = 'https://docs.google.com/spreadsheet/ccc?key=0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc#gid=1'

url2 = 'https://docs.google.com/spreadsheet/pub?key=0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc&output=html'

url3 = 'https://docs.google.com/spreadsheet/pub?key=0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc&output=html#gid=1'

exp_mapping_lines = """#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tTreatment\tDOB\tDescription
#Example mapping file for the QIIME analysis package.  These 9 samples are from a study of the effects of exercise and diet on mouse cardiac physiology (Crawford, et al, PNAS, 2009).
PC.354\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._354
PC.355\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\tControl_mouse_I.D._355
PC.356\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061126\tControl_mouse_I.D._356
PC.481\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20070314\tControl_mouse_I.D._481
PC.593\tAGCAGCACTTGT\tYATGCTGCCTCCCGTAGGAGT\tControl\t20071210\tControl_mouse_I.D._593
PC.607\tAACTGTGCGTAC\tYATGCTGCCTCCCGTAGGAGT\tFast\t20071112\tFasting_mouse_I.D._607
PC.634\tACAGAGTCGGCT\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\tFasting_mouse_I.D._634
PC.635\tACCGCAGAGTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\tFasting_mouse_I.D._635
PC.636\tACGGTGAGTGTC\tYATGCTGCCTCCCGTAGGAGT\tFast\t20080116\tFasting_mouse_I.D._636
"""


if __name__ == "__main__":
    main()
