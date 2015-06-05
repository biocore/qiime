#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Meg Pirrung", "Jesse Stombaugh", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Meg Pirrung"
__email__ = "meganap@gmail.com"

from os.path import split, splitext
import os

from biom import load_table

from qiime.util import parse_command_line_parameters, make_option
from qiime.util import load_qiime_config, create_dir, get_options_lookup
from qiime.format import (format_tep_file_lines, format_jnlp_file_lines,
                          format_te_prefs)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Makes TopiaryExplorer project file"
script_info['script_description'] = (
    "This script makes a TopiaryExplorer project file (.tep) and a jnlp file "
    "with the data location preloaded.\n\nWARNING: The jnlp file relies on an "
    "absolute path, if you move the .tep file, the generated jnlp will no "
    "longer work. However, you can still open the .tep file from your normal "
    "TopiaryExplorer install.")

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Example",
     "Create .tep file and .jnlp file:",
     "%prog -i otu_table.biom -m Fasting_Map.txt -t rep_set.tre"))

script_info['output_description'] = (
    "The result of this script is written to a .tep file and a .jnlp file, "
    "both with the name supplied by -o")

script_info['required_options'] = [
    options_lookup['otu_table_as_primary_input'],
    options_lookup['mapping_fp'],
    make_option(
        '-t',
        '--tree_fp',
        type="existing_filepath",
        help='path to tree')
]

script_info['optional_options'] = [
    options_lookup['output_dir'],
    make_option('-p', '--prefs_file_fp', type="existing_filepath",
                help='path to prefs file'),
    make_option('-w', '--web_flag', action='store_true', default=False,
                help='web codebase jnlp flag [default: %default]'),
    make_option('-u', '--url', type="string",
                help='url path for the tep file. Note: when passing this '
                'flag, it will overwrite the supplied OTU table, Mapping and '
                'Tree files.')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # get command line arguments
    otu_table_fp = opts.otu_table_fp
    mapping_fp = opts.mapping_fp
    tree_fp = opts.tree_fp
    output_dir = opts.output_dir
    output_basename = splitext(split(otu_table_fp)[1])[0]
    web_flag = opts.web_flag
    url = opts.url

    # create a generic output folder if not present
    if not output_dir:
        output_dir = 'make_tep_output/'
    create_dir(output_dir)

    # open files for parsing
    otu_table_data = load_table(otu_table_fp)

    mapping_lines = open(mapping_fp, 'U')
    tree_lines = open(tree_fp, 'U')
    prefs_dict = {}
    if opts.prefs_file_fp:
        prefs_fp = opts.prefs_file_fp
        # I don't like this eval statement
        prefs_dict = eval(open(prefs_fp, 'U').read())

    # get the tep file lines
    lines = format_tep_file_lines(otu_table_data, mapping_lines, tree_lines,
                                  prefs_dict)

    # write tep file lines if url path not given
    tep_fp = '%s/%s.tep' % (output_dir, output_basename)
    if not url:
        # write tep file
        tepfile = open(tep_fp, 'w')
        tepfile.writelines(lines)
        tepfile.close()

    # get the jnlp file lines
    lines = format_jnlp_file_lines(web_flag, url, tep_fp)

    # write jnlp file
    jnlp_fp = '%s/%s.jnlp' % (output_dir, output_basename)
    jnlpfile = open(jnlp_fp, 'w')
    jnlpfile.writelines(lines)
    jnlpfile.close()

if __name__ == "__main__":
    main()
