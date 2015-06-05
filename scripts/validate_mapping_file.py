#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from string import letters, digits

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option, create_dir
from qiime.check_id_map import check_mapping_file

options_lookup = get_options_lookup()
script_info = {}
script_info['brief_description'] = """Checks user's metadata mapping file for \
required data, valid format"""
script_info['script_description'] = """Specifically, we check that:

    - The BarcodeSequence, LinkerPrimerSequences, and ReversePrimer fields
       have valid IUPAC DNA characters, and BarcodeSequence characters
       are non-degenerate (error)
    - The SampleID, BarcodeSequence, LinkerPrimerSequence, and Description
       headers are present. (error)
    - There are not duplicate header fields (error)
    - There are not duplicate barcodes (error)
    - Barcodes are of the same length.  Suppressed when
       variable_len_barcode flag is passed (warning)
    - The headers do not contain invalid characters (alphanumeric and
       underscore only) (warning)
    - The data fields do not contain invalid characters (alphanumeric,
       underscore, space, and +-%./:,; characters) (warning)
    - SampleID fields are MIENS compliant (only alphanumeric
       and . characters). (warning)
    - There are no duplicates when the primer and variable length
       barcodes are appended (error)
    - There are no duplicates when barcodes and added demultiplex
       fields (-j option) are combined (error)
    - Data fields are not found beyond the Description column (warning)

    Details about the metadata mapping file format can be found here:
    http://www.qiime.org/documentation/file_formats.html#metadata-mapping-files

    Errors and warnings are saved to a log file.  Errors can be caused by
    problems with the headers, invalid characters in barcodes or primers, or
    by duplications in SampleIDs or barcodes.

    Warnings can arise from invalid characters and variable length barcodes that
    are not specified with the --variable_len_barcode.
    Warnings will contain a reference to the cell (row,column) that the
    warning arose from.

    In addition to the log file, a "corrected_mapping" file will be created.
    Any invalid characters will be replaced with '.' characters in
    the SampleID fields (to enforce MIENS compliance) and text in other data
    fields will be replaced with the character specified by the -c parameter,
    which is an underscore "_" by default.

    A html file will be created as well, which will show locations of
    warnings and errors, highlighted in yellow and red respectively.  If no
    errors or warnings were present the file will display a message saying
    such.  Header errors can mask other errors, so these should be corrected
    first.

    If pooled primers are used, separate with a comma.  For instance, a pooled
    set of three 27f primers (used to increase taxonomic coverage) could be
    specified in the LinkerPrimerSequence fields as such:
    AGGGTTCGATTCTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTGGTTCAG
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""", """Check the Fasting_Map.txt \
    mapping file for problems, supplying the required mapping file, and output \
    the results in the validate_mapping_file_output directory""", """%prog -m \
    Fasting_Map.txt -o validate_mapping_file_output"""))
script_info['output_description'] = """A log file, html file, and \
corrected_mapping.txt file will be written to the current output directory."""
script_info['required_options'] = [
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='Metadata mapping filepath')

]
script_info['optional_options'] = [
    make_option('-o', '--output_dir', type='new_dirpath',
                help='Required output directory for log file, corrected mapping ' +
                'file, and html file. [default: %default]', default="./"),
    make_option('-v', '--verbose',
                help='Enable printing information to standard out ' +
                '[default: %default]', default=True, action='store_false'),
    make_option('-c', '--char_replace', type='string',
                help='Changes the default character used to replace invalid ' +
                'characters found in the mapping file.  Must be a valid character (' +
                'alphanumeric, period, or underscore).' +
                '[default: %default]', default="_"),
    make_option('-b', '--not_barcoded',
                action='store_true', default=False,
                help='Use -b if barcodes are not present.  BarcodeSequence header ' +
                'still required.  [default: %default]'),
    make_option('-B', '--variable_len_barcodes',
                action='store_true', default=False,
                help='Use -B if variable length barcodes are present to suppress ' +
                'warnings about barcodes of unequal length. [default: %default]'),
    make_option('-p', '--disable_primer_check',
                action='store_true', default=False,
                help='Use -p to disable checks for primers.  LinkerPrimerSequence ' +
                'header still required. [default: %default]'),
    make_option('-j', '--added_demultiplex_field', type='string',
                help='Use -j to add a field to use in the mapping file as ' +
                'additional demultiplexing (can be used with or without barcodes).  ' +
                'All combinations ' +
                'of barcodes/primers and the these fields must be unique. The ' +
                'fields must contain values that can be parsed from the fasta labels ' +
                'such as "plate=R_2008_12_09".  In this case, "plate" would be the ' +
                'column header and "R_2008_12_09" would be the field data (minus ' +
                'quotes) in the mapping file.  To use the run prefix from the fasta ' +
                'label, such as ">FLP3FBN01ELBSX", where "FLP3FBN01" is generated ' +
                'from the run ID, use "-j run_prefix" and set the run prefix to ' +
                'be used as the data under the column header "run_prefix". ' +
                ' [default: %default]'),
    make_option('-s', '--suppress_html',
                action='store_true', default=False,
                help='Use -s to disable html file generation, can be useful for ' +
                'extremely large mapping files. [default: %default]'), ]

script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(suppress_verbose=True, **script_info)

    mapping_fp = opts.mapping_fp
    has_barcodes = not opts.not_barcoded
    variable_len_barcodes = opts.variable_len_barcodes
    output_dir = opts.output_dir + "/"
    char_replace = opts.char_replace
    verbose = opts.verbose
    disable_primer_check = opts.disable_primer_check
    added_demultiplex_field = opts.added_demultiplex_field
    suppress_html = opts.suppress_html

    # Create output directory, check path/access to mapping file
    create_dir(output_dir)

    # Test for valid replacement characters
    valid_replacement_chars = digits + letters + "_" + "."
    if char_replace not in valid_replacement_chars:
        option_parser.error('-c option requires alphanumeric, period, or ' +
                            'underscore character.')
    if len(char_replace) != 1:
        option_parser.error('-c parameter must be a single character.')

    check_mapping_file(mapping_fp, output_dir, has_barcodes, char_replace,
                       verbose, variable_len_barcodes,
                       disable_primer_check, added_demultiplex_field, suppress_html)


if __name__ == "__main__":
    main()
