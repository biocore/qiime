#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Contains functionality to interact with remote services."""

from collections import defaultdict
from csv import writer
from re import sub
from socket import gaierror
from StringIO import StringIO
from burrito.util import ApplicationNotFoundError


def raise_gdata_not_found_error(*args, **kwargs):
    raise ApplicationNotFoundError("gdata cannot be found.\nIs it installed? "
                                   "Is it in your $PYTHONPATH?\nThis is an optional QIIME "
                                   "dependency, but is required if you plan to use QIIME's remote "
                                   "mapping file features. For more information, please see "
                                   "http://qiime.org/install/install.html.")

# Load gdata if it's available. If it's not, skip it but set up to raise errors
# if the user tries to use it.
try:
    from gdata.spreadsheet import SpreadsheetsCellsFeedFromString
    from gdata.spreadsheet.service import CellQuery
    from gdata.spreadsheet.service import SpreadsheetsService
except ImportError:
    # Set functions which cannot be imported to raise_gdata_not_found_error.
    SpreadsheetsCellsFeedFromString = CellQuery = SpreadsheetsService = \
        raise_gdata_not_found_error


class GoogleSpreadsheetError(Exception):
    pass


class GoogleSpreadsheetConnectionError(Exception):
    pass


def load_google_spreadsheet(spreadsheet_key, worksheet_name=None):
    """Downloads and exports a Google Spreadsheet in TSV format.

    Returns a string containing the spreadsheet contents in TSV format (e.g.
    for writing out to a file or parsing).

    The first line is assumed to be the spreadsheet header (i.e. containing
    column names), which can optionally be followed by one or more comment
    lines (starting with '#'). Only the first cell of a comment line will be
    parsed (to keep exported spreadsheets consistent with QIIME mapping files'
    comments). The (optional) comments section is then followed by the
    spreadsheet data.

    Some of this code is based on the following websites, as well as the
    gdata.spreadsheet.text_db module:
        http://www.payne.org/index.php/Reading_Google_Spreadsheets_in_Python
        http://stackoverflow.com/a/12031835

    Arguments:
        spreadsheet_key - the key used to identify the spreadsheet (a string).
            Can either be a key or a URL containing the key
        worksheet_name - the name of the worksheet to load data from (a
            string). If not supplied, will use first worksheet in the
            spreadsheet
    """
    spreadsheet_key = _extract_spreadsheet_key_from_url(spreadsheet_key)
    gd_client = SpreadsheetsService()

    try:
        worksheets_feed = gd_client.GetWorksheetsFeed(spreadsheet_key,
                                                      visibility='public',
                                                      projection='basic')
    except gaierror:
        raise GoogleSpreadsheetConnectionError("Could not establish "
                                               "connection with server. Do "
                                               "you have an active Internet "
                                               "connection?")

    if len(worksheets_feed.entry) < 1:
        raise GoogleSpreadsheetError("The Google Spreadsheet with key '%s' "
                                     "does not have any worksheets associated "
                                     "with it." % spreadsheet_key)

    # Find worksheet that will be exported. If a name has not been provided,
    # use the first worksheet.
    worksheet = None
    if worksheet_name is not None:
        for sheet in worksheets_feed.entry:
            if sheet.title.text == worksheet_name:
                worksheet = sheet

        if worksheet is None:
            raise GoogleSpreadsheetError("The worksheet name '%s' could not "
                                         "be found in the Google Spreadsheet "
                                         "with key '%s'."
                                         % (worksheet_name, spreadsheet_key))
    else:
        # Choose the first one.
        worksheet = worksheets_feed.entry[0]

    # Extract the ID of the worksheet.
    worksheet_id = worksheet.id.text.split('/')[-1]

    # Now that we have a spreadsheet key and worksheet ID, we can read the
    # data. First get the headers (first row). We need this in order to grab
    # the rest of the actual data in the correct order (it is returned
    # unordered).
    headers = _get_spreadsheet_headers(gd_client, spreadsheet_key,
                                       worksheet_id)
    if len(headers) < 1:
        raise GoogleSpreadsheetError("Could not load spreadsheet header (it "
                                     "appears to be empty). Is your Google "
                                     "Spreadsheet with key '%s' empty?"
                                     % spreadsheet_key)

    # Loop through the rest of the rows and build up a list of data (in the
    # same row/col order found in the spreadsheet).
    spreadsheet_lines = _export_spreadsheet(gd_client, spreadsheet_key,
                                            worksheet_id, headers)

    out_lines = StringIO()
    tsv_writer = writer(out_lines, delimiter='\t', lineterminator='\n')
    tsv_writer.writerows(spreadsheet_lines)
    return out_lines.getvalue()


def _extract_spreadsheet_key_from_url(url):
    """Extracts a key from a URL in the form '...key=some_key&foo=42...

    If the URL doesn't look valid, assumes the URL is the key and returns it
    unmodified.
    """
    result = url

    if 'key=' in url:
        result = url.split('key=')[-1].split('#')[0].split('&')[0]

    return result


def _get_spreadsheet_headers(client, spreadsheet_key, worksheet_id):
    """Returns a list of headers (the first line of the spreadsheet).

    Will be in the order they appear in the spreadsheet.
    """
    headers = []

    query = CellQuery()
    query.max_row = '1'
    query.min_row = '1'
    feed = client.GetCellsFeed(spreadsheet_key, worksheet_id, query=query,
                               visibility='public', projection='values')

    # Wish python had a do-while...
    while True:
        for entry in feed.entry:
            headers.append(entry.content.text)

        # Get the next set of cells if needed.
        next_link = feed.GetNextLink()

        if next_link:
            feed = client.Get(next_link.href,
                              converter=SpreadsheetsCellsFeedFromString)
        else:
            break

    return headers


def _export_spreadsheet(client, spreadsheet_key, worksheet_id, headers):
    """Returns a list of lists containing the entire spreadsheet.

    This will include the header, any comment lines, and the spreadsheet data.
    Blank cells are represented as None. Data will only be read up to the first
    blank line that is encountered (this is a limitation of the Google
    Spreadsheet API).

    Comments are only supported after the header and before any real data is
    encountered. The lines must start with [optional whitespace] '#' and only
    the first cell is kept in that case (to avoid many empty cells after the
    comment cell, which mimics QIIME's mapping file format).

    Only cell data that falls under the supplied headers will be included.
    """
    # Convert the headers into Google's internal "cleaned" representation.
    # These will be used as lookups to pull out cell data.
    cleaned_headers = _get_cleaned_headers(headers)

    # List feed skips header and returns rows in the order they appear in the
    # spreadsheet.
    spreadsheet_lines = [headers]
    rows_feed = client.GetListFeed(spreadsheet_key, worksheet_id,
                                   visibility='public', projection='values')
    while True:
        found_data = False

        for row in rows_feed.entry:
            line = []

            # Loop through our headers and use the cleaned version to look up
            # the cell data. In certain cases (if the original header was blank
            # or only contained special characters) we will not be able to map
            # our header, so the best we can do is tell the user to change the
            # name of their header to be something simple/alphanumeric.
            for header_idx, (header, cleaned_header) in \
                    enumerate(zip(headers, cleaned_headers)):
                try:
                    cell_data = row.custom[cleaned_header].text
                except KeyError:
                    raise GoogleSpreadsheetError("Could not map header '%s' "
                                                 "to Google Spreadsheet's internal representation "
                                                 "of the header. We suggest changing the name of "
                                                 "the header in your Google Spreadsheet to be "
                                                 "alphanumeric if possible, as this will likely "
                                                 "solve the issue. Note that the name isn't "
                                                 "*required* to be alphanumeric, but it may fix "
                                                 "issues with converting to Google Spreadsheet's "
                                                 "internal format in some cases." % header)

                # Special handling of comments (if it's a comment, only keep
                # that cell to avoid several blank cells following it).
                if not found_data and header_idx == 0 and \
                   cell_data.lstrip().startswith('#'):
                    line.append(cell_data)
                    break
                else:
                    line.append(cell_data)
                    found_data = True

            spreadsheet_lines.append(line)

        # Get the next set of rows if necessary.
        next_link = rows_feed.GetNextLink()

        if next_link:
            rows_feed = client.Get(next_link.href,
                                   converter=SpreadsheetsListFeedFromString)
        else:
            break

    return spreadsheet_lines


def _get_cleaned_headers(headers):
    """Creates a list of "cleaned" headers which spreadsheets accept.

    A Google Spreadsheet converts the header names into a "cleaned" internal
    representation, which must be used to reference a cell at a particular
    header/column. They are all lower case and contain no spaces or special
    characters. If two columns have the same name after being sanitized, the
    columns further to the right have _2, _3 _4, etc. appended to them.

    If there are column names which consist of all special characters, or if
    the column header is blank, an obfuscated value will be used for a column
    name. This method does not handle blank column names or column names with
    only special characters.

    Taken from gdata.spreadsheet.text_db.ConvertStringsToColumnHeaders and
    modified to handle headers with pound signs or that start with numbers, as
    well as correctly handle duplicate cleaned headers.
    """
    cleaned_headers = []
    for header in headers:
        # Google strips special characters, whitespace, and underscores first,
        # and then strips any *leading* digits. This order is extremely
        # important!
        sanitized = sub(r'^\d+', '', sub(r'[\W_]', '', header.lower()))
        if len(sanitized) > 0:
            cleaned_headers.append(sanitized)
        else:
            raise GoogleSpreadsheetError("Encountered a header '%s' that was "
                                         "either blank or consisted only of special characters. "
                                         "Could not map the header to the internal representation "
                                         "used by the Google Spreadsheet. Please change the header "
                                         "to consist of at least one alphanumeric character."
                                         % header)

    # When the same sanitized header appears multiple times in the first row
    # of a spreadsheet, _n is appended to the name to make it unique.
    header_count = defaultdict(int)
    results = []

    for header, cleaned_header in zip(headers, cleaned_headers):
        new_header = cleaned_header

        if header_count[cleaned_header] > 0:
            # Google's numbering starts from _2, hence the +1.
            new_header = '%s_%d' % (cleaned_header,
                                    header_count[cleaned_header] + 1)

        header_count[cleaned_header] += 1
        results.append(new_header)

    return results
