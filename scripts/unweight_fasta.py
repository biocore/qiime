#!/usr/bin/env python
# File created on 20 Jun 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from skbio.parse.sequences import parse_fasta
from qiime.util import make_option
from qiime.util import parse_command_line_parameters


script_info = {}
script_info[
    'brief_description'] = "Transform fasta files with abundance weighting into unweighted"
script_info['script_description'] = '''E.g. makes 3 fasta records from a weighted input fasta file containing the following record:
>goodsample1_12_3 bc_val=20
AATGCTTGTCACATCGATGC
'''
script_info['script_usage'] = [("", '''make 3 fasta records from the following record:
>goodsample1_12_3 bc_val=20
AATGCTTGTCACATCGATGC

resulting in:
>goodsample_0
AATGCTTGTCACATCGATGC
>goodsample_1
AATGCTTGTCACATCGATGC
>goodsample_2
AATGCTTGTCACATCGATGC''', "%prog -i input.fna -o output.fna -l goodsample")]
script_info['output_description'] = "a .fasta file"
script_info['required_options'] = [
    make_option(
        '-i',
        '--input_fasta',
        type='existing_filepath',
        help='the input fasta file'),
    make_option(
        '-o',
        '--output_file',
        type='new_filepath',
        help='the output fasta filepath'),
    make_option(
        '-l',
        '--label',
        type='string',
        help='sequence label used for all records. fasta label lines will look like: >label_423'),
]
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    seqsfna_fh = open(opts.output_file, 'w')
    seq_counter = 0
    for label, seq in parse_fasta(open(opts.input_fasta, 'U')):
        seq_abundance = int(label.split()[0].split('_')[-1])
        for i in range(seq_abundance):  # don't use i, use seq_counter
            seqsfna_fh.write('>' + opts.label + '_' + str(seq_counter) + '\n')
            seqsfna_fh.write(seq + '\n')
            seq_counter += 1


if __name__ == "__main__":
    main()
