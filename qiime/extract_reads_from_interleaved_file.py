#!/usr/bin/env python

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

from cogent.parse.fastq import MinimalFastqParser
from os.path import join

from qiime.util import qiime_open


def extract_reads_from_interleaved(
        input_fp, single_file_identifier, output_dir):
    """Parses a single fasta file and creates two new files: forward and reverse, based on
    the two values (comma separated) in single_file_identifier

    input_fp: file path to input
    single_file_identifier: comma separated values to identify forward and reverse reads
    index_in_single_file: boolean defining if the indexes are also in the input_fp
    output_folder: file path to the output folder
    """
    values = single_file_identifier.split(',')
    if len(values) != 2:
        raise ValueError('You need to have two values for the headers separated by '
                         'comma but you passed "%s"' % single_file_identifier)

    forward_id = values[0]
    reverse_id = values[1]
    forward_fp = join(output_dir, "forward_reads.fastq")
    reverse_fp = join(output_dir, "reverse_reads.fastq")
    ffp = open(forward_fp, 'w')
    rfp = open(reverse_fp, 'w')

    for label, seq, qual in MinimalFastqParser(qiime_open(input_fp), strict=False):
        fastq_string = '@%s\n%s\n+\n%s\n' % (label, seq, qual)
        if forward_id in label and reverse_id not in label:
            ffp.write(fastq_string)
        elif reverse_id in label and forward_id not in label:
            rfp.write(fastq_string)
        else:
            ffp.close()
            rfp.close()
            raise ValueError("One of the input sequences doesn't have either identifier "
                             "or it has both.\nLabel: %s\nForward: %s\n Reverse: %s" %
                             (label, forward_id, reverse_id))
    ffp.close()
    rfp.close()
