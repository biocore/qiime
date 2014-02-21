#!/usr/bin/env python

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option,
                        create_dir)
from qiime.extract_reads_from_interleaved_file import extract_reads_from_interleaved

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Extract reads from an interleaved file."""
script_info[
    'script_description'] = """This script takes an interleaved file, like the ones produced by JGI, and outputs a forward and reverse fastq file with the corresponding reads in each file. """
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Extract reads from an interleaved file:""",
     """""",
     """ %prog -i $PWD/reads_to_extract.fastq -o $PWD/extracted_reads"""))
script_info[
    'output_description'] = """A new folder with two fastq files: forward_reads.fastq and reverse_reads.fastq"""
script_info['required_options'] = [
    make_option('-i', '--input_fp', type="existing_filepath",
                help='Path to input forward reads in FASTQ format.'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='Directory to store result files'),
]
script_info['optional_options'] = [
    make_option('--read_direction_identifiers', type="string",
                help='This is the string identifying the forward and reverse reads. In '
                'general, the format is 1:N:0 for read one (forward) and 2:N:0 for read '
                'two (reverse). Each value should be passed by a comma separated value. '
                '[default: %default].', default="1:N:0,2:N:0"),
]

script_info['version'] = __version__


def main():
    # parse command line parameters
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create local copy of options
    input_fp = opts.input_fp
    output_dir = opts.output_dir
    read_direction_identifiers = opts.read_direction_identifiers

    create_dir(output_dir, fail_on_exist=False)

    extract_reads_from_interleaved(
        input_fp,
        read_direction_identifiers,
        output_dir)


if __name__ == "__main__":
    main()
