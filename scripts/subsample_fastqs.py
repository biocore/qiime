#!/usr/bin/env python

__author__ = "Jose Antonio Navas Molina"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jose Antonio Navas Molina"
__email__ = "josenavasmolina@gmail.com"
__status__ = "Development"

from qiime.util import parse_command_line_parameters, make_option, \
    subsample_fastqs

script_info={}
script_info['brief_description']="""Randomly subsample sequences from two fastq\
 files"""
script_info['script_description']="""Subsample the reads and barcode files of \
a given file at the same time."""
script_info['script_usage']=[]
script_info['script_usage'].append(
    ("""Example:""",
    """Subsample seqs.fastq and seqs_barcodes.fastq to approximately 5%""",
    """%prog -i $PWD/seqs.fastq -b $PWD/seqs_barcodes.fastq -p 0.05 """ +\
    """-o $PWD/subsampled_seqs.fastq -u $PWD/barcodes_subsampled.fastq"""))
script_info['output_description']=""""""
script_info['required_options']=[
    make_option('-i', '--input_fp', type='existing_filepath',
        help='Path to the fastq file.'),
    make_option('-b', '--barcodes_fp', type='existing_filepath',
        help='Path to the barcodes fastq file.'),
    make_option('-p','--percent_subsample',action='store',type='float',
        help='Specify the percentage of sequences to subsample')
]
script_info['optional_options']=[\
    make_option('-o', '--output_fp', type='new_filepath',
        help='Path to the output fastq file'),
    make_option('-u', '--barcodes_output_fp', type='new_filepath',
        help='Path to the output barcode fastq file')
]
script_info['version'] = __version__

if __name__ == '__main__':
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_fp = opts.input_fp
    barcodes_fp = opts.barcodes_fp
    percent_subsample = opts.percent_subsample
    output_fp = opts.output_fp
    barcodes_output_fp = opts.barcodes_output_fp

    if percent_subsample > 1 or percent_subsample <= 0:
        raise ValueError, 'percent_subsample must be in range 0-1'

    if not output_fp:
        in_file_basename, in_file_ext = splitext(split(input_fp)[1])
        output_fp = '%s_subsample_%3.2f%s' % (in_file_basename,
            percent_subsample, in_file_ext)

    if not barcodes_output_fp:
        in_file_basename, in_file_ext = splitext(split(barcodes_fp)[1])
        barcodes_output_fp = '%s_subsample_%3.2f%s' % (in_file_basename,
            percent_subsample, in_file_ext)

    subsample_fastqs(input_fp, output_fp, barcodes_fp, barcodes_output_fp,
        percent_subsample)