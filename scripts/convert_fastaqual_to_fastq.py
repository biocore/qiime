#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Development"

from qiime.util import make_option, create_dir,\
 parse_command_line_parameters, get_options_lookup
from qiime.convert_fastaqual_to_fastq import convert_fastq

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description'] = "From a FASTA file and a matching QUAL file, generates a minimal FASTQ file."
script_info['script_description'] = "From a FASTA file and a mathcing QUAL file, generates a minimal FASTQ file. A minimal FASTQ file omits the redundtant sequence label on the quality scores; the quality scores for a sequence are assumed to follow immediately after the sequence with which they are associated. The output FASTQ file will be generated in the specified output directory with the same name as the input FASTA file, suffixed with '.fastq'"

script_info['script_usage'] = []
script_info['script_usage'].append(("Example:",\
        "Using the input files seqs.fna and seqs.qual, generate seqs.fastq in the fastq_files directory:", "python fastaQualToFastq_script.py -f seqs.fna -q seqs.qual -o fastq_files/"))

script_info['output_description'] = """Outputs a minimal FASTQ file, which omits the redundant sequence label on the quality scores."""
script_info['required_options'] = [
            make_option('-f', '--fasta_fp',
                        type='existing_filepath',
                                help='Input FASTA file.'),

            make_option('-q', '--qual_fp', type='existing_filepath',
                help='Input QUAL file.')
            ]

script_info['optional_options'] = [
            make_option('-o', '--output_dir',
                        type = 'new_path',
                        help = 'Output directory. Will be created if does not \
exist. [default: %default]', default="."),

            make_option('-a', '--ascii_increment',
                type = int,
                help = 'The number to add to the quality score to get the \
ASCII character. [default: %default]', default = 33),

            make_option('-F', '--full_fasta_headers',
                action = 'store_true',
                help = 'Include full FASTA headers in FASTQ file \
(as opposed to merely the sequence label). [default: %default]', \
                default = False),

            make_option('-i', '--full_fastq',
                action = 'store_true',
                help = 'Include identifiers on quality lines in the FASTQ file \
(those beginning with a "+" [default=%default]', default = False),

            make_option('-m', '--multiple_output_files',
                action = 'store_true',
                help = 'Create multiple FASTQ files, one for each sample. \
[default=%default]', default = False)
            ]
        
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    fasta_fp = opts.fasta_fp
    qual_fp = opts.qual_fp
    output_dir = opts.output_dir
    full_fastq = opts.full_fastq
    ascii_increment = opts.ascii_increment
    full_fasta_headers = opts.full_fasta_headers
    multiple_output_files = opts.multiple_output_files
    create_dir(output_dir)
    convert_fastq(fasta_fp, qual_fp, output_dir, multiple_output_files, 
            ascii_increment, full_fastq, full_fasta_headers)

if __name__ == "__main__":
    main()
