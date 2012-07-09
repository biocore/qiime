#!/usr/bin/env python

__author__ = "Adam Robbins-Pianka, Abhisaar Yadav"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Adam Robbins-Pianka, Abhisaar Yadav"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Adam Robbins-Pianka"
__email__ = "adam.robbinspianka@colorado.edu"
__status__ = "Development"

# Reviewed by William Walters

from os import path

from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser

from qiime.parse import parse_qual_score

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

    fasta_file = open(fasta_file_path,'U')
    qual_file = open(qual_file_path,'U')
    
    if not multiple_output_files:
		output_file_path = path.join(output_directory, \
		path.splitext(path.split(fasta_file_path)[1])[0] + '.fastq')
		try:
			fastq_file = open(output_file_path,'w')
		except IOError:
			qual_file.close()
			fasta_file.close()
			raise IOError,("Could not open FASTQ file for writing: " \
					+ output_file_path + '\n')
    output_files = {}
    
    # Reading in entire files to get SampleIDs.  Also can't assume that the
    # qual and fasta files are in the same order, although they should be.
    # This should probably be a separate function.
    
    # Read in entire FASTA file
    fasta = [x for x in MinimalFastaParser(fasta_file)]
    fasta_file.close()

    # Read in entire QUAL file
    qual = parse_qual_score(qual_file)
    qual_file.close()
    
    sample_ids = []
    for header, sequence in fasta:
        label = header.split()[0]
        sample_id = label.split('_')[0]
        sample_ids.append(sample_id)
    
    sample_id_counter = {}
    sample_id_counts = {}
    for k in sample_ids:
        sample_id_counts[k] = sample_ids.count(k)
        sample_id_counter[k] = 0
    
    for header, sequence in fasta:
        label = header.split()[0]
        sample_id = label.split('_')[0]
        sample_id_counter[sample_id]+=1
        
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
                except IOError:
                    raise IOError,("Could not open FASTQ file for writing: " \
                            + output_file_path + '\n')
            fastq_file = output_files[sample_id]
            
        try:
            quality_scores = qual[label]
        except KeyError:
            raise KeyError,("No entry in QUAL file for label: %s\n" % label)

        if len(quality_scores) != len(sequence):
            raise KeyError,("Number of quality scores "+\
            "(%d) does not match number of positions (%d) for label: %s" %\
             (len(quality_scores), len(sequence), label))
                
        #Writing to FASTQ file
        fastq_file.write('@' + fastq_sequence_header + '\n')
        fastq_file.write(sequence + '\n')
        fastq_file.write('+' + fastq_quality_header + '\n')
        for qual_score in qual[label]:
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
    
    fastq_file = open(fasta_file_path,'U')  
    
    if not multiple_output_files:
		output_fasta = path.join(output_directory, \
			path.splitext(path.split(fastq_file_path)[1])[0] + '.fna')
		output_qual = path.join(output_directory, \
			path.splitext(path.split(fastq_file_path)[1])[0] + '.qual')
				
		try:
			fasta_o = open(output_fasta,'w')
		except IOError:
			fastq_file.close()
			raise IOError,("Could not open FASTA file for writing: %s" %\
			 output_fasta + '\n')
		
		try:
			qual_o = open(output_qual,'w')
		except IOError:
			fastq_file.close()
			raise IOError,("Could not open QUAL file for writing: %s" %\
			 output_qual + '\n')
    fasta_output = {}
    qual_output = {}
        
    #Read FASTQ file    
    fastq = [y for y in MinimalFastqParser(fastq_file, strict=False)]
    fastq_file.close() 
    
    sample_ids = []
    for header, sequence, qual in fastq:
        label = header.split()[0]
        sample_id = label.split('_')[0]
        sample_ids.append(sample_id)
    
    sample_id_counter = {}
    sample_id_counts = {}
    for k in sample_ids:
        sample_id_counts[k] = sample_ids.count(k)
        sample_id_counter[k] = 0
        
    for header, sequence, qual in fastq:
        label = header.split()[0]
        sample_id = label.split('_')[0]
        sample_id_counter[sample_id]+=1
        
        if full_fasta_headers: label = header
        
        if multiple_output_files:
            if sample_id not in fasta_output:
                output_fasta = path.join(output_directory, \
                    path.splitext(path.split(fastq_file_path)[1])[0] + \
                     '_' + sample_id + '.fna')
                try:
                    fasta_output[sample_id] = open(output_fasta,'w')
                except IOError:
                    raise IOError,("Could not open output FASTA file: %s" %\
                     output_fasta + '\n')
                    
            if sample_id not in qual_output:
                output_qual = path.join(output_directory, \
                 path.splitext(path.split(fastq_file_path)[1])[0] +'_'+ \
                 sample_id +'.qual')
                try:
                    qual_output[sample_id] = open(output_qual,'w')
                except IOError:
                    fastq_file.close()
                    raise IOError,("Could not open QUAL file for writing: %s" %\
                     output_qual + '\n')
                     
            fasta_o = fasta_output[sample_id]
            qual_o = qual_output[sample_id]
        
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
