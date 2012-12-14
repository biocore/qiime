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

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)
from qiime.remote import load_google_spreadsheet_mapping_file

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Downloads and saves a remote mapping file"
script_info['script_description'] = """
This script exports, downloads, and saves a mapping file that is stored
remotely. Currently, the only type of remote mapping file that is supported is
a Google Spreadsheet, though other methods of remote storage may be supported
in the future.

For more information and examples pertaining to this script and remote mapping
files in general, please refer to the accompanying tutorial, which can be found
at http://qiime.org/tutorials/remote_mapping_files.html.
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("", "", "%prog"))

script_info['output_description'] = """
The script outputs a single file, which is the metadata mapping file obtained
from the remote location (in QIIME-compatible format).
"""

script_info['required_options'] = [
    make_option('-k', '--spreadsheet_key', type='string',
        help='the spreadsheet key that will be used to identify the Google '
        'Spreadsheet to load. This is the part of the Google Spreadsheet URL '
        'that comes after \'key=\'. You may instead provide the entire URL '
        'and the key will be extracted from it. Either this option or '),
    options_lookup['output_fp']
]
script_info['optional_options'] = [
    make_option('-w', '--worksheet_name', type='string',
        help='the name of the worksheet in the Google Spreadsheet that '
        'contains the mapping file. If the worksheet name contains spaces, '
        'please include quotes around the name. [default: the first worksheet '
        'in the Google Spreadsheet will be used]', default=None)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    worksheet_name = opts.worksheet_name
    if worksheet_name is not None:
        worksheet_name = worksheet_name.strip('"').strip("'")

    results = load_google_spreadsheet_mapping_file(opts.spreadsheet_key,
                                                   worksheet_name)

    with open(opts.output_fp, 'w') as out_f:
        out_f.write(results)


if __name__ == "__main__":
    main()
