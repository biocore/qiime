#!/usr/bin/env python
# file: join_paired_ends.py

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Mike Robeson"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Mike Robeson"
__email__ = "robesonms@ornl.gov"
__status__ = "Development"

from cogent.parse.fastq import MinimalFastqParser
from cogent.app.fastq_join import join_paired_end_reads_fastqjoin
from cogent.app.seqprep import join_paired_end_reads_seqprep
import os
import gzip

join_method_constructors = {}
join_method_names = {'fastq-join':join_paired_end_reads_fastqjoin,
                     'SeqPrep':join_paired_end_reads_seqprep}


def write_synced_barcodes_fastq(joined_fp, index_fp):
    """Writes new index file based on surviving assembled paired-ends.
       -joined_fp : file path to paired-end assembled fastq file
       -index_fp : file path to index / barcode reads fastq file

       This function iterates through the joined reads file and index file. 
       Only those index-reads within the file at index_fp, that have headers
       matching those within the joined-pairs at joined_fp, are written 
       to file. 

     WARNING: Assumes reads are in the same order in both files,
              except for cases in which the corresponding
              read in the joined_fp file is missing (i.e. pairs 
              failed to assemble).

    """

    j_path, j_ext = os.path.splitext(joined_fp)
    i_path, i_ext = os.path.splitext(index_fp)

    # check if joined pairs file is gzipped
    # if not, read normally
    if j_ext == '.gz':
        jh = gzip.open(joined_fp,'r')
    else:
        jh = open(joined_fp,'U')
    
    # check if index / barcode file is gzipped
    # if not, read normally
    if i_ext == '.gz':
        ih = gzip.open(index_fp,'r')
    else:
        ih = open(index_fp,'U')
    
    # base new index file name on joined paired-end file name:
    filtered_bc_outfile_path = j_path + '_barcodes.fastq'
    fbc_fh = open(filtered_bc_outfile_path, 'w')


    # Set up iterators
    index_fastq_iter = MinimalFastqParser(ih, strict=False)
    joined_fastq_iter = MinimalFastqParser(jh, strict=False) 
    # Write barcodes / index reads that we observed within
    # the joined paired-ends. Warn if index and joined data
    # are not in order.
    for joined_label,joined_seq,joined_qual in joined_fastq_iter:
        index_label,index_seq,index_qual = index_fastq_iter.next()
        while joined_label != index_label:
            try:
                index_label,index_seq,index_qual = index_fastq_iter.next()
            except StopIteration:
                raise StopIteration, "\n\nReached end of index-reads file"+\
                 " before iterating through joined paired-end-reads file!"+\
                 " Except for missing paired-end reads that did not survive"+\
                 " assembly, your index and paired-end reads files must be in"+\
                 " the same order! Last ID processed was:\n"+ \
                 " \'%s\'\n" %(joined_label)
        else:
            fastq_string = '@%s\n%s\n+\n%s\n'\
                            %(index_label,index_seq,index_qual)
            fbc_fh.write(fastq_string)
    
    ih.close()
    jh.close()
    fbc_fh.close()

    return filtered_bc_outfile_path

        
        
