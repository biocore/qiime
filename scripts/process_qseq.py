#!/usr/bin/env python
# File created on 30 Mar 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import make_option
from glob import glob
from os.path import split
from skbio.util import create_dir
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.format import illumina_data_to_fastq

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Given a directory of per-swath qseq files,\
 this script generates a single fastq per lane."
script_info['script_description'] = ""
script_info['script_usage'] = [
    ("", "Generate fastq files from all lanes of read 1 data in the current\
 directory.", "process_qseq.py -i ./ -o ./fastq/ -r 1"),
    ("", "Generate fastq files from all lanes of read 2 data in the current\
 directory, truncating the sequences after the first 12 bases.",
     "process_qseq.py -i ./ -o ./fastq/ -r 2 -b 12")]
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_dir', type='existing_dirpath',
                help='the input directory'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output directory'),
    make_option('-r', '--read', help='the read number to consider', type='int')
]
script_info['optional_options'] = [
    make_option('-l', '--lanes', type='string',
                help='the lane numbers to consider, comma-separated [defaut: %default]',
                default='1,2,3,4,5,6,7,8'),
    make_option('-b', '--bases', type='int',
                help='the number of bases to include (useful for slicing a barcode)\
 [defaut: all]',
                default=None),
    make_option('--ignore_pass_filter', action='store_true', default=False,
                help='ignore the illumina pass filter [default:%default; reads with 0 in ' +
                ' pass filter field are discarded]')
]
script_info['version'] = __version__


def iter_split_lines(lines):
    """ """
    for line in lines:
        yield tuple(line.split('\t'))


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    input_dir = opts.input_dir
    output_dir = opts.output_dir
    create_dir(output_dir)
    lanes = opts.lanes.split(',')
    bases = opts.bases
    read = opts.read
    ignore_pass_filter = opts.ignore_pass_filter

    for lane in lanes:
        read1_fps = sorted(glob('%s/s_%s_%d_*qseq.txt' % (input_dir,
                                                          lane.replace(
                                                              ',',
                                                              ''),
                                                          read)))
        # sort so results will be consistent across different runs (important
        # so amplicon and barcodes read headers will match)
        output_fp = '%s/s_%s_%s_sequences.fastq' % (output_dir, lane, read)
        output_f = open(output_fp, 'w')
        for read1_fp in read1_fps:
            for record in iter_split_lines(open(read1_fp, 'U')):
                fastq_s, pass_filter = illumina_data_to_fastq(record,
                                                              number_of_bases=bases)
                if ignore_pass_filter or pass_filter != 0:
                    output_f.write('%s\n' % fastq_s)
        output_f.close()

if __name__ == "__main__":
    main()
