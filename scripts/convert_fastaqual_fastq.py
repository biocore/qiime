#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka, Abhisaar Yadav"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka, Abhisaar Yadav",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

# Reviewed by William Walters

from qiime.util import make_option, create_dir,\
    parse_command_line_parameters, get_options_lookup
from qiime.convert_fastaqual_fastq import convert_fastaqual_fastq

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "From a FASTA file and a matching QUAL file,\
 generates a FASTQ file. From FASTQ file generates FASTA file and \
 matching QUAL file."
script_info['script_description'] = "From a FASTA file and a matching QUAL \
file, generates a FASTQ file. A minimal FASTQ file omits the redundant \
sequence label on the quality scores; the quality scores for a sequence are \
assumed to follow immediately after the sequence with which they are \
associated. The output FASTQ file will be generated in the specified output \
directory with the same name as the input FASTA file, suffixed with '.fastq'. \
A FASTQ file will be split into FASTA and QUAL files, and generated in the \
designated output directory."

script_info['script_usage'] = []
script_info['script_usage'].append(("Example:",
                                    "Using the input files seqs.fna and seqs.qual, generate seqs.fastq in \
the fastq_files directory:", "%prog -f seqs.fna -q seqs.qual -o fastq_files/"))
script_info['script_usage'].append(("Example:",
                                    "Using input seqs.fastq generate fasta and qual files in fastaqual \
directory:", "%prog -c fastq_to_fastaqual \
-f seqs.fastq -o fastaqual"))
script_info['output_description'] = """Outputs a complete or minimal FASTQ \
file, which omits the redundant sequence label on the quality scores, or splits\
 FASTQ file into matching FASTA/QUAL files."""

script_info['required_options'] = [
    make_option('-f', '--fasta_file_path',
                type='existing_filepath',
                help='Input FASTA or FASTQ file.')]

script_info['optional_options'] = [

    make_option('-q', '--qual_file_path', type='existing_filepath',
                help='Required input QUAL file if converting to FASTQ.',
                default=None),

    make_option('-o', '--output_dir',
                type='new_path',
                help='Output directory. Will be created if does not ' +
                'exist. [default: %default]', default="."),

    make_option('-c', '--conversion_type',
                type='choice', choices=['fastaqual_to_fastq', 'fastq_to_fastaqual'],
                help='type of conversion: fastaqual_to_fastq or ' +
                'fastq_to_fastaqual [default: %default]', default=
                "fastaqual_to_fastq"),

    make_option('-a', '--ascii_increment',
                type='int',
                help='The number to add (subtract if coverting from FASTQ) ' +
                'to the quality score to get the ASCII character (or numeric ' +
                'quality score). [default: %default]', default=33),

    make_option('-F', '--full_fasta_headers',
                action='store_true',
                help='Include full FASTA headers in output file(s) (as ' +
                'opposed to merely the sequence label). [default: %default]',
                default=False),

    make_option('-b', '--full_fastq',
                action='store_true',
                help='Include identifiers on quality lines in the FASTQ ' +
                'file (those beginning with a "+"). Irrelevant when ' +
                'converting from FASTQ. [default=%default]', default=False),

    make_option('-m', '--multiple_output_files',
                action='store_true',
                help='Create multiple FASTQ files, one for each sample, or ' +
                'create multiple matching FASTA/QUAL for each sample. ' +
                '[default=%default]', default=False)]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    fasta_file_path = opts.fasta_file_path
    qual_file_path = opts.qual_file_path
    output_dir = opts.output_dir
    full_fastq = opts.full_fastq
    ascii_increment = opts.ascii_increment
    full_fasta_headers = opts.full_fasta_headers
    multiple_output_files = opts.multiple_output_files
    conversion_type = opts.conversion_type

    create_dir(output_dir)

    convert_fastaqual_fastq(fasta_file_path, qual_file_path, conversion_type,
                            output_dir, multiple_output_files, ascii_increment, full_fastq,
                            full_fasta_headers)

if __name__ == "__main__":
    main()
