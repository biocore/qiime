#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.5.0"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Release"


from sys import stderr
from os import path
from cogent.parse.fasta import MinimalFastaParser
from qiime.parse import parse_qual_score

# TODO: unit testing!
def convert_fastq(fasta_file_path, qual_file_path, output_directory='.',
        multiple_output_files=False, ascii_increment=33,
        full_fastq=False, full_fasta_headers=False):
    '''Takes a FASTA file and it's corresponding QUAL file, and returns
a FASTQ file. Returns True on success, False on failure.'''
    # Make sure the input files are accessible
    try:
        qual_file = open(qual_file_path)
    except IOError, e:
        stderr.write("Could not open QUAL file: " + qual_file_path + '\n')
        return False

    try:
        fasta_file = open(fasta_file_path)
    except IOError, e:
        qual_file.close()
        stderr.write("Could not open FASTA file: " + fasta_file_path + '\n')
        return False

    if not multiple_output_files:
        output_file_path = path.join(output_directory, \
                path.splitext(path.split(fasta_file_path)[1])[0] + '.fastq')
        try:
            fastq_file = open(output_file_path,'w')
        except IOError, e:
            qual_file.close()
            fasta_file.close()
            stderr.write("Could not open FASTQ file for writing: " \
                    + output_file_path + '\n')
            return False
    output_files = {}

    # Read in entire FASTA file
    fasta = [x for x in MinimalFastaParser(fasta_file)]
    fasta_file.close()

    # Read in entire QUAL file
    qual = parse_qual_score(qual_file)
    qual_file.close()

    for header, sequence in fasta:
        label = header.split()[0]
        sample_id = label.split('_')[0]

        fastq_sequence_header = label
        if full_fasta_headers: fastq_sequence_header = header
        fastq_quality_header = ''
        if full_fastq: fastq_quality_header = fastq_sequence_header

        if multiple_output_files:
            if sample_id not in output_files:
                output_file_path = path.join(output_directory, \
                        path.splitext(path.split(fasta_file_path)[1])[0] + \
                        '_' + sample_id + '.fastq')
                try:
                    output_files[sample_id] = open(output_file_path,'w')
                except IOError, e:
                    stderr.write("Could not open FASTQ file for writing: " \
                            + output_file_path + '\n')
                    return False
            fastq_file = output_files[sample_id]

        try:
            quality_scores = qual[label]
        except KeyError, e:
            stderr.write("No entry in QUAL file for label: %s\n" % label)
            fastq_file.close()
            if multiple_output_files:
                for k in output_files.keys():
                    if not output_files[k].closed: output_files[k].close()
            return False

        if len(quality_scores) != len(sequence):
            stderr.write("Number of quality scores (%d) does not match number \
of positions (%d) for label: %s" % (len(quality_scores), len(sequence), label))
            fastq_file.close()
            if multiple_output_files:
                for k in output_files.keys():
                    if not output_files[k].closed: output_files[k].close()
            return False

        fastq_file.write('@' + fastq_sequence_header + '\n')
        fastq_file.write(sequence + '\n')
        fastq_file.write('+' + fastq_quality_header + '\n')
        for qual_score in qual[label]:
            # increment the qual score by the asciiIncrement (default 33),
            # and print the corresponding character, which represents that 
            # position's quality.
            qual_score += ascii_increment
            if qual_score < 32 or qual_score > 126:
                raise ValueError, "Cannot convert quality score to ASCII code \
between 32 and 126: " + str(qual_score - ascii_increment) + " using \
ascii_increment = " + str(ascii_increment)
            fastq_file.write(chr(qual_score))
        fastq_file.write('\n')
    fastq_file.close()
    if multiple_output_files:
        for k in output_files.keys():
            if not output_files[k].closed: output_files[k].close()
    return True
