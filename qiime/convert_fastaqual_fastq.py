#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka, Abhisaar Yadav"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka, Abhisaar Yadav"]
__license__ = "GPL"
__version__ = "1.6.0"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Release"

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
    """ Checks files/permissions, calls appropriate conversion function 
    
    fasta_file_path:  filepath of input FASTA or FASTQ file.
    qual_file_path:  filepath of input QUAL file (needed for making FASTQ files)
    conversion_type:  Either fastqual_to_fastq or fastq_to_fastqual.
    output_directory:  Directory to output converted files.
    multiple_output_files:  Make one file per SampleID.
    ascii_increment:  Conversion value for fastq ascii character to numeric
     quality score.
    full_fastq:  Write labels to both sequence and quality score lines.
    full_fasta_headers:  Retain all data on fasta label, instead of breaking at
     first whitespace."""
    
    try:
        fasta_file = open(fasta_file_path,'U')
        fasta_file.close()
    except IOError:
        raise IOError,("Could not open FASTA/FASTQ file: " +\
        fasta_file_path + '\n')
    fasta_file.close()
    
    if conversion_type == 'fastaqual_to_fastq':
        if qual_file_path == None:
            raise ValueError,("Must specify a QUAL filepath when converting "+\
             "to fastq format.")
        try:
            qual_file = open(qual_file_path,'U')
            qual_file.close()
        except IOError:
            raise IOError,("Could not open QUAL file: " + qual_file_path + '\n')
        convert_fastq(fasta_file_path, qual_file_path, output_directory,
        multiple_output_files, ascii_increment,
        full_fastq, full_fasta_headers);

    elif conversion_type == 'fastq_to_fastaqual':
        convert_fastaqual(fasta_file_path, output_directory,
        multiple_output_files, ascii_increment,
        full_fastq, full_fasta_headers);

    else:
        raise ValueError,('conversion_type must be fastaqual_to_fastq '+ \
         'or fastq_to_fastaqual.')

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
    
    
    output_files = {}
    
    fasta_file = open(fasta_file_path,'U')
    qual_file = open(qual_file_path,'U')
    
    
    
    # Need to open file the first time as "w", thereafter open as "a"
    sample_ids_written = {}
    
    for fasta_data, qual_data in izip(MinimalFastaParser(fasta_file),
         MinimalQualParser(qual_file)):
        
        qual_header = qual_data[0]
        fasta_header = fasta_data[0] 
        label = fasta_header.split()[0]
        sample_id = label.split('_')[0]
        sequence = fasta_data[1]
        qual = qual_data[1]
        try: quality_scores = qual_data[1]
        except KeyError:
            raise KeyError,("No entry in QUAL file for label: %s\n" % \
            label)
            
        if qual_header != label:
            raise KeyError,("Fasta(%s) and qual(%s) headers don't match" %\
            (label, qual_header))
            
        if len(qual) != len(sequence):
            raise KeyError,("Number of quality scores "+\
            "(%d) does not match number of positions (%d) for label: %s" %\
             (len(qual), len(sequence), label))

        
            
        if not multiple_output_files:
            output_file_path = path.join(output_directory, \
            path.splitext(path.split(fasta_file_path)[1])[0] + '.fastq')
            if output_file_path in sample_ids_written.keys():
                sample_ids_written[output_file_path] = True
            else:
                sample_ids_written[output_file_path] = False
            try:
                # Create new file if first time writing, else append
                if sample_ids_written[output_file_path]:
                    fastq_file = open(output_file_path, 'a')
                else:
                    fastq_file = open(output_file_path, 'w')
            except IOError:
                qual_file.close()
                fasta_file.close()
                raise IOError,("Could not open FASTQ file for writing: " \
                        + output_file_path + '\n')
            output_files[sample_id] = output_file_path
                
        if multiple_output_files:
            if sample_id not in output_files:
                output_file_path = path.join(output_directory, \
                        path.splitext(path.split(fasta_file_path)[1])[0] + \
                        '_' + sample_id + '.fastq')
                if output_file_path in sample_ids_written.keys():
                    sample_ids_written[output_file_path] = True
                else:
                    sample_ids_written[output_file_path] = False
                try:
                    # Create new file if first time writing, else append
                    if sample_ids_written[output_file_path]:
                        output_files[sample_id] = open(output_file_path, 'a')
                    else:
                        output_files[sample_id] = open(output_file_path, 'w')
                    
                except IOError:
                    raise IOError,("Could not open FASTQ file for writing: " \
                            + output_file_path + '\n')
                output_files[sample_id] = output_file_path

        fastq_file = open(output_files[sample_id], 'a')

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
        qual_scores = list(qual)
        for qual_score in qual_scores:
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
        if multiple_output_files:
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
    
    fastq_file_path = fasta_file_path
    
    fasta_output = {}
    qual_output = {}
    fastq_file = open(fasta_file_path,'U')
    
    # Need to open file the first time as "w", thereafter open as "a"
    sample_ids_written = {}
    
    for fastq_data in izip(MinimalFastqParser(fastq_file, strict=False)):
        sequence = fastq_data[0][1]
        qual = fastq_data[0][2]
        header = fastq_data[0][0]
        label = header.split()[0]
        sample_id = label.split('_')[0]
        
        if len(sequence) != len(qual):
            raise KeyError,("Number of quality scores "+\
            "(%d) does not match number of positions (%d) for label: %s" %\
             (len(qual), len(sequence), label))

    
        if not multiple_output_files:
            output_fasta = path.join(output_directory, \
                path.splitext(path.split(fastq_file_path)[1])[0] + '.fna')
            output_qual = path.join(output_directory, \
                path.splitext(path.split(fastq_file_path)[1])[0] + '.qual')
                
            if output_fasta in sample_ids_written.keys():
                sample_ids_written[output_fasta] = True
            else:
                sample_ids_written[output_fasta] = False
            try:
                # Create new file if first time writing, else append
                if sample_ids_written[output_fasta]:
                    fasta_o = open(output_fasta,'a')
                    qual_o = open(output_qual,'a')
                else:
                    fasta_o = open(output_fasta,'w')
                    qual_o = open(output_qual,'w')
            except IOError:
                raise IOError,("Could not open output FASTA or QUAL files, "+\
                 "please check file permissions.")

            fasta_output[sample_id] = output_fasta
            qual_output[sample_id] = output_qual
            
        if multiple_output_files:
            if sample_id not in fasta_output:
                output_fasta = path.join(output_directory, \
                    path.splitext(path.split(fastq_file_path)[1])[0] + \
                     '_' + sample_id + '.fna')
                     
                if output_fasta in sample_ids_written.keys():
                    sample_ids_written[output_fasta] = True
                else:
                    sample_ids_written[output_fasta] = False
                    
                try:
                    if sample_ids_written[output_fasta]:
                        fasta_output[sample_id] = open(output_fasta, 'a')
                    else:
                        fasta_output[sample_id] = open(output_fasta, 'w')
                except IOError:
                    raise IOError,("Could not open output FASTA file: %s" %\
                     output_fasta + '\n')
                fasta_output[sample_id] = output_fasta
                
            if sample_id not in qual_output:
                output_qual = path.join(output_directory, \
                 path.splitext(path.split(fastq_file_path)[1])[0] +'_'+ \
                 sample_id +'.qual')
                try:
                    if sample_ids_written[output_fasta]:
                        qual_output[sample_id] = open(output_qual, 'a')
                    else:
                        qual_output[sample_id] = open(output_qual, 'w')
                    
                    #qual_output[sample_id] = open(output_qual,'a')
                except IOError:
                    fastq_file.close()
                    raise IOError,("Could not open QUAL file for writing: %s" %\
                     output_qual + '\n')
                qual_output[sample_id] = output_qual
                
        if full_fasta_headers: label = header
        
        fasta_o = open(fasta_output[sample_id], 'a')
        qual_o = open(qual_output[sample_id], 'a')
        
        #write Fasta file   
        fasta_o.write('>' + label + '\n') 
        fasta_o.write(sequence + '\n')

        #convert quality scores
        qual_chars = list(qual)
        qual_scores = []
        for qual_char in qual_chars:
            if (ord(qual_char) - ascii_increment) < -0: 
                raise ValueError,("Output qual scores are negative values. "+ \
                 "Use different ascii_increment value than %s" %\
                 str(ascii_increment))
            else:
                qual_scores.append(ord(qual_char) - ascii_increment)
        
        #write QUAL file        
        score_numbers = []
        for i, qual_score in enumerate(qual_scores):
            score_numbers.append(i)
        qual_o.write('>' + label +'\n')        
        for i, qual_score in enumerate(qual_scores):
            if i % 60 == 0 and i != 0:
                qual_o.write('\n')
            qual_o.write(str(qual_score))
            if (i+1) % 60 != 0 and i != max(score_numbers):
                qual_o.write(' ')
        qual_o.write('\n')
        if multiple_output_files:
            fasta_o.close()
            qual_o.close()
