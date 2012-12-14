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

"""Contains functionality to interact with remote services."""

from csv import writer
from socket import gaierror
from StringIO import StringIO
from gdata.spreadsheet import SpreadsheetsCellsFeedFromString
from gdata.spreadsheet.service import CellQuery
from gdata.spreadsheet.service import SpreadsheetsService

class RemoteMappingFileError(Exception):
    pass

# TODO test comments, empty lines/cells, quoted strings
def load_google_spreadsheet_mapping_file(spreadsheet_key, worksheet_name=None):
    """Loads a mapping file contained in a Google Spreadsheet.

    Returns a string containing the mapping file contents in QIIME-compatible
    format (e.g. for writing out to a file or parsing using
    qiime.parse.parse_mapping_file).

    Some of this code is based on the following websites, as well as the
    gdata.spreadsheet.text_db module:
        http://www.payne.org/index.php/Reading_Google_Spreadsheets_in_Python
        http://stackoverflow.com/a/12031835
    """
    gd_client = SpreadsheetsService()

    try:
        worksheets_feed = gd_client.GetWorksheetsFeed(spreadsheet_key,
                                                      visibility='public',
                                                      projection='basic')
    except gaierror:
        raise RemoteMappingFileError("Could not establish connection with "
                                     "server. Do you have an active Internet "
                                     "connection?")

    if len(worksheets_feed.entry) < 1:
        raise RemoteMappingFileError("The Google Spreadsheet with key '%s' "
                                     "does not have any worksheets associated "
                                     "with it." % spreadsheet_key)

    # Find worksheet that will be used as the mapping file. If a name has not
    # been provided, use the first worksheet.
    worksheet = None
    if worksheet_name is not None:
        for sheet in worksheets_feed.entry:
            if sheet.title.text == worksheet_name:
                worksheet = sheet

        if worksheet is None:
            raise RemoteMappingFileError("The worksheet name '%s' could not "
                                         "be found in the Google Spreadsheet "
                                         "with key '%s'."
                                         % (worksheet_name, spreadsheet_key))
    else:
        # Choose the first one.
        worksheet = worksheets_feed.entry[0]

    # Extract the ID of the worksheet.
    worksheet_id = worksheet.id.text.split('/')[-1]

    # Now that we have a spreadsheet key and worksheet ID, we can read the
    # mapping file data. First get the mapping file headers (first row). We
    # need this in order to grab the rest of the actual mapping file data in
    # the correct order (it is returned unordered).
    query = CellQuery()
    query.max_row = '1'
    query.min_row = '1'
    feed = gd_client.GetCellsFeed(spreadsheet_key, worksheet_id, query=query,
                                  visibility='public', projection='values')

    # Wish python had a do-while...
    headers = []
    while True:
        for entry in feed.entry:
            headers.append(entry.content.text)

        # Get the next set of cells if needed.
        next_link = feed.GetNextLink()

        if next_link:
            feed = gd_client.Get(next_link.href,
                                 converter=SpreadsheetsCellsFeedFromString)
        else:
            break

    if len(headers) < 1:
        raise RemoteMappingFileError("Could not load mapping file header (it "
                                     "appears to be empty). Is your Google "
                                     "Spreadsheet with key '%s' empty?"
                                     % spreadsheet_key)

    # Convert the actual headers into Google's internal "cleaned"
    # representation.
    cleaned_headers = _convert_strings_to_column_headers(headers)

    # Loop through the rest of the rows and build up a list of data (in the
    # same row/col order found in the original mapping file).
    rows_feed = gd_client.GetListFeed(spreadsheet_key, worksheet_id,
                                      visibility='public', projection='values')
    mapping_lines = [headers]
    while True:
        for row in rows_feed.entry:
            try:
                mapping_lines.append([row.custom[cleaned_header].text
                    for header, cleaned_header in zip(headers, cleaned_headers)])
            except KeyError:
                raise RemoteMappingFileError("Could not map header '%s' to Google "
                                             "Spreadsheet's internal "
                                             "representation of the header. We "
                                             "suggest changing the name of the "
                                             "header in your Google Spreadsheet "
                                             "to be alphanumeric if possible, as "
                                             "this will likely solve the issue." %
                                             header)

        # Get the next set of rows if necessary.
        next_link = rows_feed.GetNextLink()

        if next_link:
            rows_feed = gd_client.Get(next_link.href,
                                      converter=SpreadsheetsListFeedFromString)
        else:
            break

    out_lines = StringIO()
    tsv_writer = writer(out_lines, delimiter='\t', lineterminator='\n')
    tsv_writer.writerows(mapping_lines)
    return out_lines.getvalue()

def _extract_spreadsheet_key_from_url(url):
    """Extracts a key from a URL in the form '...key=some_key#foo=42...
    
    If the URL doesn't look valid, assumes the URL is the key and returns it
    unmodified.
    """
    result = url

    if 'docs.google.com' in url:
        result = url.split('key=')[-1].split('#')[0]

    return result

def _convert_strings_to_column_headers(proposed_headers):
  """Converts a list of strings to column names which spreadsheets accepts.

  When setting values in a record, the keys which represent column names must
  fit certain rules. They are all lower case, contain no spaces or special
  characters. If two columns have the same name after being sanitized, the 
  columns further to the right have _2, _3 _4, etc. appended to them.

  If there are column names which consist of all special characters, or if
  the column header is blank, an obfuscated value will be used for a column
  name. This method does not handle blank column names or column names with
  only special characters.

  Taken from gdata.spreadsheet.text_db.ConvertStringsToColumnHeaders and
  modified to handle headers with pound signs.
  """
  headers = []
  for input_string in proposed_headers:
    # Probably a more efficient way to do this. Perhaps regex.
    sanitized = input_string.lower().replace('_', '').replace(
        ':', '').replace(' ', '').replace('#', '')
    # When the same sanitized header appears multiple times in the first row
    # of a spreadsheet, _n is appended to the name to make it unique.
    header_count = headers.count(sanitized)
    if header_count > 0:
      headers.append('%s_%i' % (sanitized, header_count+1))
    else:
      headers.append(sanitized)
  return headers
