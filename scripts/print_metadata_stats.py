#!/usr/bin/env python
# File created on 19 Dec 2012
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Daniel McDonald", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from collections import defaultdict
from qiime.parse import parse_mapping_file
from qiime.util import parse_command_line_parameters, make_option
from qiime.sort import natsort
from sys import stdout

script_info = {}
script_info[
    'brief_description'] = "Count the number of samples associated to a category value"
script_info[
    'script_description'] = """Sum up the number of samples with each category value and print this information."""
script_info[
    'script_usage'] = [("Example:", "Count the number of samples associated with Treatment", """%prog -m $PWD/mapping.txt -c Treatment"""),
                       ("Example writting the output to a file", "Count the number of samples associated with Treatment and save them to a file called stats.txt", """%prog -m mapping.txt -c Treatment -o stats.txt""")]
script_info[
    'output_description'] = """Two columns, the first being the category value and the second being the count. Output is to standard out. If there are unspecified values, the output category is identified as ***UNSPECIFIED***"""
script_info['required_options'] = [
    make_option(
        '-m',
        '--mapping_file',
        type="existing_filepath",
        help='the input metadata file'),
    make_option(
        '-c',
        '--category',
        type='string',
        help='the category to examine')
]
script_info['optional_options'] = [
    make_option('-o', '--output_fp', type="new_filepath",
                help="path where output will be written [default: print to screen]",
                default=None)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    output_fp = opts.output_fp

    map_data, header, comments = parse_mapping_file(opts.mapping_file)

    if opts.category not in header:
        option_parser.error(
            "%s doesn't appear to exist in the mapping file!" %
            opts.category)

    # use stdout or the user supplied file path
    if output_fp:
        fd = open(output_fp, 'w')
    else:
        fd = stdout

    result = defaultdict(int)
    cat_idx = header.index(opts.category)
    for samp in map_data:
        result[samp[cat_idx]] += 1

    for cat_val in natsort(result):
        if not cat_val:
            fd.write("***UNSPECIFIED***\t%d\n" % result[cat_val])
        else:
            fd.write("%s\t%d\n" % (cat_val, result[cat_val]))

    fd.close()

if __name__ == "__main__":
    main()
