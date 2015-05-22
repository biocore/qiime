#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka, Abhisaar Yadav"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka, Abhisaar Yadav"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"

# Reviewed by William Walters

from os import path
from itertools import izip
from collections import defaultdict

from qiime.parse import QiimeParseError, MinimalQualParser
from skbio.parse.sequences import parse_fasta
from skbio.parse.sequences import parse_fastq


def convert_fastaqual_fastq(fasta_file_path, qual_file_path,
                            conversion_type='fastaqual_to_fastq', output_directory='.',
                            multiple_output_files=False, ascii_increment=33,
                            full_fastq=False, full_fasta_headers=False):
    """Calls appropriate conversion function, depending on direction.

    fasta_file_path:  filepath of input FASTA or FASTQ file.
    qual_file_path:  filepath of input QUAL file (needed for making FASTQ files)
    conversion_type:  Either fastqual_to_fastq or fastq_to_fastqual.
    output_directory:  Directory to output converted files.
    multiple_output_files:  Make one file per SampleID.
    ascii_increment:  Conversion value for fastq ascii character to numeric
     quality score.
    full_fastq:  Write labels to both sequence and quality score lines.
    full_fasta_headers:  Retain all data on fasta label, instead of breaking at
     first whitespace.
     """

    if conversion_type == 'fastaqual_to_fastq':
        convert_fastq(fasta_file_path, qual_file_path, output_directory,
                      multiple_output_files, ascii_increment,
                      full_fastq, full_fasta_headers)

    elif conversion_type == 'fastq_to_fastaqual':
        convert_fastaqual(fasta_file_path, output_directory,
                          multiple_output_files, ascii_increment,
                          full_fastq, full_fasta_headers)

    else:
        raise ValueError('conversion_type must be fastaqual_to_fastq '
                         'or fastq_to_fastaqual.')


def get_filename_with_new_ext(original_file_path, new_ext, output_directory):
    """Returns the original file name, but with a different extension

    E.g.
    get_filname_with_new_ext('/Users/shared/test.fasta', '.fastq', '.')
    returns
    'test.fastq'
    """
    return path.join(output_directory,
                     path.splitext(path.split(original_file_path)[1])[0] + new_ext)


def convert_fastq(fasta_file_path, qual_file_path, output_directory='.',
                  multiple_output_files=False, ascii_increment=33,
                  full_fastq=False, full_fasta_headers=False,
                  per_file_buffer_size=100000):
    '''Takes a FASTA and QUAL file, generates FASTQ file(s)

    fasta_file_path:  filepath of input FASTA file.
    qual_file_path:  filepath of input QUAL file (needed for making FASTQ files)
    output_directory:  Directory to output converted files.
    multiple_output_files:  Make one file per SampleID.
    ascii_increment:  Conversion value for fastq ascii character to numeric
     quality score.
    full_fastq:  Write labels to both sequence and quality score lines.
    full_fasta_headers:  Retain all data on fasta label, instead of breaking at
     first whitespace.'''

    fasta_file = open(fasta_file_path, 'U')
    qual_file = open(qual_file_path, 'U')

    # if we're not using multiple output files, we can open the one (and only)
    # output file right now
    if not multiple_output_files:
        output_file_path = get_filename_with_new_ext(fasta_file_path,
                                                     '.fastq',
                                                     output_directory)

        fastq_file = open(output_file_path, 'w')
    else:
        fastq_lookup = defaultdict(str)

    # iterate through the FASTA and QUAL files entry by entry (assume the
    # entries are synchronized)
    for fasta_data, qual_data in izip(parse_fasta(fasta_file),
                                      MinimalQualParser(qual_file)):

        qual_header = qual_data[0]
        fasta_header = fasta_data[0]

        label = fasta_header.split()[0]
        sample_id = label.split('_')[0]

        sequence = fasta_data[1]
        qual = qual_data[1]

        # check whether the entries are actually (at least nominally) synch'd
        if qual_header != label:
            raise KeyError(("QUAL header (%s) does not match "
                            "FASTA header (%s)") % (qual_header, label))

        if len(sequence) != len(qual):
            raise KeyError(("Sequence length does not match QUAL length for "
                            "label (%s)") % label)

        if multiple_output_files:
            output_file_path = get_filename_with_new_ext(fasta_file_path,
                                                         '_' + sample_id +
                                                         '.fastq',
                                                         output_directory)

            # when we use multiple output files, we close each file after each
            # sequence is written to avoid using up all the file handles, so
            # we must open the file each time in append mode
            # fastq_file = open(output_file_path, 'a')

        if full_fasta_headers:
            fastq_sequence_header = fasta_header
        else:
            fastq_sequence_header = label

        if full_fastq:
            fastq_quality_header = fastq_sequence_header
        else:
            fastq_quality_header = ''

        # Writing to FASTQ file
        record = '@%s\n%s\n+%s\n' % (fastq_sequence_header,
                                     sequence,
                                     fastq_quality_header)

        if multiple_output_files:
            fastq_lookup[output_file_path] += record
        else:
            fastq_file.write(record)

        for qual_score in qual:
            # increment the qual score by the asciiIncrement (default 33),
            # and print the corresponding character, which represents that
            # position's quality.
            qual_score += ascii_increment
            if qual_score < 32 or qual_score > 126:
                raise ValueError("Cannot convert quality score to ASCII code" +
                                 " between 32 and 126: " + str(qual_score - ascii_increment) +
                                 "using ascii_increment = " + str(ascii_increment))

            if multiple_output_files:
                fastq_lookup[output_file_path] += chr(qual_score)
            else:
                fastq_file.write(chr(qual_score))

        if multiple_output_files:
            fastq_lookup[output_file_path] += '\n'
        else:
            fastq_file.write('\n')

        if multiple_output_files:
            if len(fastq_lookup[output_file_path]) >= per_file_buffer_size:
                fastq_file = open(output_file_path, 'a')
                fastq_file.write(fastq_lookup[output_file_path])
                fastq_lookup[output_file_path] = ''
                fastq_file.close()

    # write last seqs to output files, or close the output file if thre is only
    # one
    if multiple_output_files:
        for output_file_path, records in fastq_lookup.iteritems():
            if records:
                fastq_file = open(output_file_path, 'a')
                fastq_file.write(records)
                fastq_file.close()
    else:
        fastq_file.close()


def convert_fastaqual(fasta_file_path, output_directory='.',
                      multiple_output_files=False, ascii_increment=33,
                      full_fastq=False, full_fasta_headers=False,
                      per_file_buffer_size=100000):
    '''Takes a FASTQfile, generates FASTA and QUAL file(s)

    fasta_file_path:  filepath of input FASTQ file.
    output_directory:  Directory to output converted files.
    multiple_output_files:  Make one file per SampleID.
    ascii_increment:  Conversion value for fastq ascii character to numeric
     quality score.
    full_fastq:  Write labels to both sequence and quality score lines.
    full_fasta_headers:  Retain all data on fasta label, instead of breaking at
     first whitespace.'''

    # rename this to avoid confusion...
    fastq_fp = fasta_file_path

    # if we are NOT using multiple output files, then open our two (and only)
    # output files here
    if not multiple_output_files:
        fasta_out_fp = get_filename_with_new_ext(fastq_fp,
                                                 '.fna',
                                                 output_directory)
        qual_out_fp = get_filename_with_new_ext(fastq_fp,
                                                '.qual',
                                                output_directory)

        fasta_out_f = open(fasta_out_fp, 'w')
        qual_out_f = open(qual_out_fp, 'w')

    else:
        fasta_out_lookup = defaultdict(str)
        qual_out_lookup = defaultdict(str)

    fpo = ascii_increment
    for header, sequence, qual in parse_fastq(open(fastq_fp, 'U'),
                                              strict=False,
                                              phred_offset=fpo):
        label = header.split()[0]
        sample_id = label.split('_')[0]

        if multiple_output_files:
            fasta_out_fp = get_filename_with_new_ext(fastq_fp,
                                                     '_' + sample_id + '.fna',
                                                     output_directory)

            qual_out_fp = get_filename_with_new_ext(fastq_fp,
                                                    '_' + sample_id + '.qual',
                                                    output_directory)

        if full_fasta_headers:
            label = header

        if (qual < 0).any():
            raise ValueError("Output qual scores are negative values. "
                             "Use different ascii_increment value than %s" %
                             str(ascii_increment))

        # write QUAL file, 60 qual scores per line
        qual_record = [">%s\n" % label]
        for i in range(0, len(qual), 60):
            qual_record.append(' '.join([str(q) for q in qual[i:i + 60]]))
            qual_record.append('\n')
        qual_record = ''.join(qual_record)

        if multiple_output_files:
            qual_out_lookup[qual_out_fp] += qual_record
        else:
            qual_out_f.write(qual_record)

        # write FASTA file
        fasta_record = '>%s\n%s\n' % (label, sequence)
        if multiple_output_files:
            fasta_out_lookup[fasta_out_fp] += fasta_record
        else:
            fasta_out_f.write(fasta_record)

        # if we're writing multiple output files, we must close after each
        # sequeunce write to avoid potentiallyusing up all the OS's filehandles
        if multiple_output_files:
            if fasta_out_lookup[fasta_out_fp] >= per_file_buffer_size:
                fasta_f = open(fasta_out_fp, 'a')
                fasta_f.write(fasta_out_lookup[fasta_out_fp])
                fasta_f.close()
                fasta_out_lookup[fasta_out_fp] = ''

                qual_f = open(qual_out_fp, 'a')
                qual_f.write(qual_out_lookup[qual_out_fp])
                qual_f.close()
                qual_out_lookup[qual_out_fp] = ''

    # if we have one output file, close it now
    if multiple_output_files:
        for fasta_out_fp, records in fasta_out_lookup.iteritems():
            if records:
                fasta_f = open(fasta_out_fp, 'a')
                fasta_f.write(records)
                fasta_f.close()

        for qual_out_fp, records in qual_out_lookup.iteritems():
            if records:
                qual_f = open(qual_out_fp, 'a')
                qual_f.write(records)
                qual_f.close()
    else:
        fasta_out_f.close()
        qual_out_f.close()
