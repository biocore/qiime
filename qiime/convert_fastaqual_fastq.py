#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka, Abhisaar Yadav"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka, Abhisaar Yadav"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Development"

# Reviewed by William Walters

from os import path


from qiime.parse import QiimeParseError, MinimalQualParser
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from itertools import izip
from qiime.parse import parse_qual_score
from time import time

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
        raise ValueError,('conversion_type must be fastaqual_to_fastq '
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
        full_fastq=False, full_fasta_headers=False):
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
    
    
    fasta_file = open(fasta_file_path,'U')
    qual_file = open(qual_file_path,'U')
    
    # if we're not using multiple output files, we can open the one (and only)
    # output file right now
    if not multiple_output_files:
        output_file_path = get_filename_with_new_ext(fasta_file_path,
                                                     '.fastq',
                                                     output_directory)

        fastq_file = open(output_file_path, 'w')

    # iterate through the FASTA and QUAL files entry by entry (assume the
    # entries are synchronized)
    for fasta_data, qual_data in izip(MinimalFastaParser(fasta_file),
         MinimalQualParser(qual_file)):
        
        qual_header = qual_data[0]
        fasta_header = fasta_data[0] 

        label = fasta_header.split()[0]
        sample_id = label.split('_')[0]

        sequence = fasta_data[1]
        qual = qual_data[1]

        # check whether the entries are actually (at least nominally) synch'd
        if qual_header != label:
            raise KeyError, ("QUAL header (%s) does not match "
                             "FASTA header (%s)") % (qual_header, label)

        if len(sequence) != len(qual):
            raise KeyError, ("Sequence length does not match QUAL length for "
                             "label (%s)") % label

        if multiple_output_files:
            output_file_path = get_filename_with_new_ext(fasta_file_path,
                                                 '_' + sample_id + '.fastq',
                                                 output_directory)

            # when we use multiple output files, we close each file after each
            # sequence is written to avoid using up all the file handles, so
            # we must open the file each time in append mode
            fastq_file = open(output_file_path, 'a')

        if full_fasta_headers:
            fastq_sequence_header = fasta_header
        else:
            fastq_sequence_header = label

        if full_fastq:
            fastq_quality_header = fastq_sequence_header
        else:
            fastq_quality_header = ''

        #Writing to FASTQ file
        fastq_file.write('@' + fastq_sequence_header + '\n')
        fastq_file.write(sequence + '\n')
        fastq_file.write('+' + fastq_quality_header + '\n')

        for qual_score in qual:
            # increment the qual score by the asciiIncrement (default 33),
            # and print the corresponding character, which represents that 
            # position's quality.
            qual_score += ascii_increment
            if qual_score < 32 or qual_score > 126:
                raise ValueError,("Cannot convert quality score to ASCII code"+\
                 " between 32 and 126: " + str(qual_score - ascii_increment) +\
                 "using ascii_increment = " + str(ascii_increment))
            fastq_file.write(chr(qual_score))

        fastq_file.write('\n')

        # Must close the output file here to avoid potentially using up all
        # the OS's filehandles
        if multiple_output_files:
            fastq_file.close()

    # if we have only one output file, close it here
    if not multiple_output_files:
        fastq_file.close()
        
def convert_fastaqual(fasta_file_path, output_directory='.',
        multiple_output_files=False, ascii_increment=33,
        full_fastq=False, full_fasta_headers=False):
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
        fasta_out_lookup = {}
        qual_out_lookup = {}

    for header, sequence, qual in MinimalFastqParser(open(fastq_fp, 'U'),
                                                     strict=False):
        label = header.split()[0]
        sample_id = label.split('_')[0]

        if multiple_output_files:
            fasta_out_fp = get_filename_with_new_ext(fastq_fp,
                                     '_' + sample_id + '.fna',
                                     output_directory)

            qual_out_fp = get_filename_with_new_ext(fastq_fp,
                                     '_' + sample_id + '.qual',
                                     output_directory)

            fasta_out_f = open(fasta_out_fp, 'a')
            qual_out_f = open(qual_out_fp, 'a')

        if full_fasta_headers:
            label = header

        #convert quality scores
        qual_scores = []
        for qual_char in qual:
            if (ord(qual_char) - ascii_increment) < 0: 
                raise ValueError,("Output qual scores are negative values. "
                 "Use different ascii_increment value than %s" %
                 str(ascii_increment))
            else:
                qual_scores.append(str(ord(qual_char) - ascii_increment))

        #write QUAL file, 60 qual scores per line
        qual_out_f.write('>' + label +'\n')
        for i in range(0, len(qual_scores), 60):
            qual_out_f.write(' '.join(qual_scores[i:i+60]) + '\n')

        #write FASTA file
        fasta_out_f.write('>' + label + '\n') 
        fasta_out_f.write(sequence + '\n')

        # if we're writing multiple output files, we must close after each
        # sequeunce write to avoid potentiallyusing up all the OS's filehandles
        if multiple_output_files:
            fasta_out_f.close()
            qual_out_f.close()

    # if we have one output file, close it now
    if not multiple_output_files:
        fasta_out_f.close()
        qual_out_f.close()
