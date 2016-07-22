#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from os.path import isdir

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option, create_dir
from qiime.add_qiime_labels import add_qiime_labels


options_lookup = get_options_lookup()
script_info = {}
script_info[
    'brief_description'] = """Takes a directory, a metadata mapping file, and a column name that contains the fasta file names that SampleIDs are associated with, combines all files that have valid fasta extensions into a single fasta file, with valid QIIME fasta labels."""
script_info['script_description'] = """A metadata mapping file with SampleIDs
and fasta file names (just the file name itself, not the full or relative
filepath) is used to generate a combined fasta file with valid
QIIME labels based upon the SampleIDs specified in the mapping file.

See: http://qiime.org/documentation/file_formats.html#metadata-mapping-files
for details about the metadata file format.

Example mapping file:
#SampleID	BarcodeSequence	LinkerPrimerSequence	InputFileName	Description
Sample.1	AAAACCCCGGGG	CTACATAATCGGRATT	seqs1.fna	sample.1
Sample.2	TTTTGGGGAAAA	CTACATAATCGGRATT	seqs2.fna	sample.2

This script is to handle situations where fasta data comes already
demultiplexed into a one fasta file per sample basis.  Only alters
the fasta label to add a QIIME compatible label at the beginning.

Example:
With the metadata mapping file above, and an specified directory containing the
files seqs1.fna and seqs2.fna, the first line from the seqs1.fna file might
look like this:
>FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_
AACAGATTAGACCAGATTAAGCCGAGATTTACCCGA

and in the output combined fasta file would be written like this
>Sample.1_0 FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_
AACAGATTAGACCAGATTAAGCCGAGATTTACCCGA

No changes are made to the sequences.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """Specify fasta_dir as the input directory of fasta files, use the metadata mapping file example_mapping.txt, with the metadata fasta file name column specified as InputFileName, start enumerating with 1000000, and output the data to the directory combined_fasta""",
     """%prog -i fasta_dir -m example_mapping.txt -c InputFileName -n 1000000 -o combined_fasta"""))
script_info[
    'output_description'] = """A combined_seqs.fasta file will be created in the output directory, with the sequences assigned to the SampleID given in the metadata mapping file."""
script_info['required_options'] = [
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='SampleID to fasta file name mapping file filepath'),
    make_option('-i', '--fasta_dir', type='existing_dirpath',
                help='Directory of fasta files to combine and label.'),
    make_option('-c', '--filename_column', type=str,
                help='Specify column used in metadata mapping file for ' +
                'fasta file names.')

]
script_info['optional_options'] = [
    make_option('-o', '--output_dir', type='new_dirpath',
                help='Required output directory for log file and corrected mapping ' +
                'file, log file, and html file. [default: %default]', default="."),
    make_option('-n', '--count_start',
                help='Specify the number to start enumerating sequence labels with. ' +
                '[default: %default]', default=0, type="int")
]

script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    mapping_fp = opts.mapping_fp
    fasta_dir = opts.fasta_dir
    output_dir = opts.output_dir
    count_start = int(opts.count_start)
    filename_column = opts.filename_column

    # Check input filepaths
    try:
        test_mapping_f = open(mapping_fp, "U")
    except IOError:
        raise IOError("Cannot open mapping filepath " +
                      "%s, please check filepath and permissions." % mapping_fp)

    if not isdir(fasta_dir):
        raise IOError("Specified fasta dir " +
                      "%s, does not exist" % fasta_dir)

    # Create output directory, check path/access to mapping file
    create_dir(output_dir)

    add_qiime_labels(open(mapping_fp, "U"), fasta_dir, filename_column,
                     output_dir, count_start)


if __name__ == "__main__":
    main()
