#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import parse_command_line_parameters, make_option, MetadataMap

script_info = {}
script_info['brief_description'] = """Merge mapping files"""
script_info[
    'script_description'] = (
        "This script provides a convenient interface for merging mapping "
        "files which contain data on different samples.")
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Example:",
     "Merge two mapping files into a new mapping file (merged_mapping.txt). "
     "In cases where a mapping field is not provided for some samples, add "
     "the value 'Data not collected'.",
     "%prog -m map_controls.txt,map_fasting.txt -o merged_mapping.txt -n "
     "'Data not collected'"))
script_info[
    'output_description'] = (
        "The result of this script is a merged mapping file (tab-delimited).")
script_info['required_options'] = [
    make_option('-m', '--mapping_fps', type='existing_filepaths',
                help='the input mapping files in a comma-separated list'),
    make_option('-o', '--output_fp', type='new_filepath',
                help='the output mapping file to write'),
]

script_info['optional_options'] = [
    make_option('-n', '--no_data_value',
                help='value to represent missing data (i.e., when all '
                'fields are not defined in all mapping files) [default: '
                '%default]', default='no_data'),
    make_option('--case_insensitive', action='store_true', default=False,
                help='if present the headers will be merged case insensitivly '
                'and transformed to upper case [default: %default]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose
    output_fp = opts.output_fp
    mapping_files = [open(fp, 'U') for fp in opts.mapping_fps]
    no_data_value = opts.no_data_value
    ci = opts.case_insensitive

    mapping_data = MetadataMap.mergeMappingFiles(mapping_files,
                                                 no_data_value=no_data_value,
                                                 case_insensitive=ci)

    with open(output_fp, 'w') as f:
        f.write(str(mapping_data))

if __name__ == "__main__":
    main()
