#!/usr/bin/env python
# File created February 29, 2012
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from qiime.util import make_option

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    create_dir
from qiime.truncate_reverse_primer import truncate_reverse_primer


options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Takes a demultiplexed fasta file, finds a\
 specified reverse primer sequence, and truncates this primer and subsequent\
 sequences following the reverse primer."""
script_info['script_description'] = """Takes input mapping file and fasta\
 sequences which have already have been demultiplexed (via split_libraries.py,\
 denoise_wrapper.py, ampliconnoise.py, etc.) with fasta labels that are in\
 QIIME format, i.e., SampleID_#. This script will use the SampleID and a\
 mapping file with a ReversePrimer column to find the reverse primer by local\
 alignment and remove this and any subsequent sequence in a filtered output\
 fasta file."""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""",
                                    """Find, truncate reverse primers from the fasta file seqs.fna, with the\
 SampleIDs and reverse primers specified in Mapping_File_Rev_Primer.txt, writes\
 output fasta file to the reverse_primer_removed directory:""",
                                    """%prog -f seqs.fna -m Mapping_File_Rev_Primer.txt\
 -o reverse_primer_removed/"""))
script_info['output_description'] = """Truncated version of the input fasta file\
 (based on input name with 'seqs_rev_primer_truncated' appended) will be\
 generated in the output directory, along with a .log file."""
script_info['required_options'] = [
    make_option('-f', '--fasta_fp',
                type='existing_filepath',
                help='Fasta file.  Needs to have fasta labels in proper ' +
                'demultiplexed format.'),

    make_option('-m', '--mapping_fp',
                type='existing_filepath',
                help="Mapping filepath.  ReversePrimer field required.  Reverse " +
                "primers need to be in 5'->3' orientation.")
]

script_info['optional_options'] = [
    make_option('-o', '--output_dir',
                type='new_path',
                help='Output directory.  Will be created if does not exist.  ' +
                '[default: %default]', default="."),
    make_option('-z', '--truncate_option', type='choice',
                choices=['truncate_only', 'truncate_remove'],
                help='Truncation option.  The default option, "truncate_only" will ' +
                'try to find the reverse primer to truncate, and if not found, ' +
                'will write the sequence unchanged.  If set to "truncate_remove", ' +
                'sequences where the reverse primer is not found will not be ' +
                'written. [default: %default]', default='truncate_only'),
    make_option('-M', '--primer_mismatches', type='int',
                help='Number of mismatches allowed in the reverse primer. ' +
                '[default: %default]', default=2)
]

script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    fasta_fp = opts.fasta_fp
    mapping_fp = opts.mapping_fp
    output_dir = opts.output_dir
    truncate_option = opts.truncate_option
    primer_mismatches = int(opts.primer_mismatches)

    create_dir(output_dir)

    if truncate_option not in ['truncate_only', 'truncate_remove']:
        raise ValueError('-z option must be either truncate_only or ' +
                         'truncate_remove')

    try:
        fasta_f = open(fasta_fp, "U")
        fasta_f.close()
    except IOError:
        raise IOError("Unable to open fasta file, please check path/" +
                      "permissions.")
    try:
        mapping_f = open(fasta_fp, "U")
        mapping_f.close()
    except IOError:
        raise IOError("Unable to open mapping file, please check path/" +
                      "permissions.")

    truncate_reverse_primer(fasta_fp, mapping_fp, output_dir, truncate_option,
                            primer_mismatches)


if __name__ == "__main__":
    main()
