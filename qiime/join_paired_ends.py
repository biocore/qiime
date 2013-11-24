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
    """Writes new index file based on surviving assembled paired-ends."""

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
    # the joined paired-ends:
    for joined_label,joined_seq,joined_qual in joined_fastq_iter:
        for index_label,index_seq,index_qual in index_fastq_iter:
            if joined_label == index_label:
                fastq_string = '@%s\n%s\n+\n%s\n'\
                                %(index_label,index_seq,index_qual)
                fbc_fh.write(fastq_string)
                break # stop looping through 2nd for-loop when match is found

    ih.close()
    jh.close()
    fbc_fh.close()



        
        
