#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the remote.py module."""

from socket import gaierror
from cogent.util.unit_test import TestCase, main
from qiime.remote import (_get_cleaned_headers,
                          _get_spreadsheet_headers,
                          _export_spreadsheet,
                          _extract_spreadsheet_key_from_url,
                          load_google_spreadsheet,
                          raise_gdata_not_found_error,
                          GoogleSpreadsheetConnectionError,
                          GoogleSpreadsheetError)

try:
    from gdata.spreadsheet.service import SpreadsheetsService
except ImportError:
    SpreadsheetsService = raise_gdata_not_found_error

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

    def getClient(self):
        """Returns a spreadsheet client if there's an Internet connection.

        Returns None if no connection.
        """
        client = SpreadsheetsService()

        try:
            client.GetWorksheetsFeed(self.spreadsheet_key, visibility='public',
                                     projection='basic')
        except gaierror:
            client = None

        return client

    def test_load_google_spreadsheet(self):
        """Test retrieves a Google Spreadsheet. Will fail if no Internet connection."""
        # Test without naming a worksheet.
        obs = load_google_spreadsheet(self.spreadsheet_key,
                                      worksheet_name=None)
        self.assertEqual(obs, self.exp_mapping_lines)

        # Test with naming a worksheet.
        obs = load_google_spreadsheet(self.spreadsheet_key,
                                      worksheet_name='Fasting_Map')
        self.assertEqual(obs, self.exp_mapping_lines)

    def test_load_google_spreadsheet_invalid_input(self):
        """Test correctly raises errors on various bad inputs. Will fail if no Internet connection."""
        # Bad worksheet name.
        self.assertRaises(GoogleSpreadsheetError, load_google_spreadsheet,
                         self.spreadsheet_key, worksheet_name='foo')

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
        """Test reading the header line from a spreadsheet. Will fail if no Internet connection."""
        client = self.getClient()
        if client:
            exp = ['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence',
                   'Treatment', 'DOB', 'Description']
            obs = _get_spreadsheet_headers(client, self.spreadsheet_key,
                                           self.worksheet_id)
            self.assertEqual(obs, exp)
        else:
            raise GoogleSpreadsheetConnectionError("Cannot execute test "
                    "without an active Internet connection.")

    def test_export_spreadsheet(self):
        """Test exporting spreadsheet as TSV. Will fail if no Internet connection."""
        client = self.getClient()
        if client:
            exp = [['#SampleID', 'DOB'],
                ['#Example mapping file for the QIIME analysis package.  '
                 'These 9 samples are from a study of the effects of exercise '
                 'and diet on mouse cardiac physiology (Crawford, et al, '
                 'PNAS, 2009).'], ['PC.354', '20061218'],
                ['PC.355', '20061218'], ['PC.356', '20061126'],
                ['PC.481', '20070314'], ['PC.593', '20071210'],
                ['PC.607', '20071112'], ['PC.634', '20080116'],
                ['PC.635', '20080116'], ['PC.636', '20080116']]
            obs = _export_spreadsheet(client, self.spreadsheet_key,
                                      self.worksheet_id, ['#SampleID', 'DOB'])
            self.assertEqual(obs, exp)
        else:
            raise GoogleSpreadsheetConnectionError("Cannot execute test "
                    "without an active Internet connection.")

    def test_export_spreadsheet_invalid_input(self):
        """Test exporting spreadsheet with bad input raises errors. Will fail if no Internet connection."""
        client = self.getClient()
        if client:
            # Nonexisting header.
            self.assertRaises(GoogleSpreadsheetError, _export_spreadsheet,
                    client, self.spreadsheet_key, self.worksheet_id,
                    ['#SampleID', 'Foo'])
        else:
            raise GoogleSpreadsheetConnectionError("Cannot execute test "
                    "without an active Internet connection.")

    def test_get_cleaned_headers(self):
        """Test correctly converts headers to Google's representation."""
        # Some duplicates.
        exp = ['foo', 'foo_2', 'foo_3', 'foo_4', 'fooo', 'foo_5', 'foo_6',
               'foo_7', 'foo_8', 'foo_9', 'f2oo456', 'foo_10']
        obs = _get_cleaned_headers(
                ['foo', 'Foo', 'FOO', 'F_oO', 'F:Oo_o', '123foo', '#Foo',
                 '123foo', ' 123Foo', 'f O\tO#', ' f2\too456', '456 foo'])
        self.assertEqual(obs, exp)

        # All unique.
        exp = ['foo', 'bar']
        obs = _get_cleaned_headers(['Fo#o', 'bar'])
        self.assertEqual(obs, exp)

        # Header consisting of only special characters and header that is
        # blank.
        self.assertRaises(GoogleSpreadsheetError, _get_cleaned_headers,
                          ['Foo', '___', 'BAR'])
        self.assertRaises(GoogleSpreadsheetError, _get_cleaned_headers,
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
