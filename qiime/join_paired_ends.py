#!/usr/bin/env python

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Mike Robeson"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Mike Robeson"
__email__ = "robesonms@ornl.gov"

from skbio.parse.sequences import parse_fastq
from skbio.format.sequences import format_fastq_record
from bfillings.fastq_join import FastqJoin, join_paired_end_reads_fastqjoin
from bfillings.seqprep import SeqPrep, join_paired_end_reads_seqprep
from qiime.util import qiime_open
import os
import gzip

join_method_constructors = {}
join_method_names = {'fastq-join': join_paired_end_reads_fastqjoin,
                     'SeqPrep': join_paired_end_reads_seqprep}


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

    # open files (handles normal / gzipped data)
    jh = qiime_open(joined_fp)
    ih = qiime_open(index_fp)

    # base new index file name on joined paired-end file name:
    j_path, ext = os.path.splitext(joined_fp)
    filtered_bc_outfile_path = j_path + '_barcodes.fastq'
    fbc_fh = open(filtered_bc_outfile_path, 'w')

    # Set up iterators
    index_fastq_iter = parse_fastq(ih, strict=False, enforce_qual_range=False)
    joined_fastq_iter = parse_fastq(jh, strict=False, enforce_qual_range=False)
    # Write barcodes / index reads that we observed within
    # the joined paired-ends. Warn if index and joined data
    # are not in order.
    for joined_label, joined_seq, joined_qual in joined_fastq_iter:
        index_label, index_seq, index_qual = index_fastq_iter.next()
        while joined_label != index_label:
            try:
                index_label, index_seq, index_qual = index_fastq_iter.next()
            except StopIteration:
                raise StopIteration("\n\nReached end of index-reads file" +
                                    " before iterating through joined paired-end-reads file!" +
                                    " Except for missing paired-end reads that did not survive" +
                                    " assembly, your index and paired-end reads files must be in" +
                                    " the same order! Also, check that the index-reads and" +
                                    " paired-end reads have identical headers. The last joined" +
                                    " paired-end ID processed was:\n\'%s\'\n" % (joined_label))
        else:
            fastq_string = format_fastq_record(index_label, index_seq, index_qual)
            fbc_fh.write(fastq_string)

    ih.close()
    jh.close()
    fbc_fh.close()

    return filtered_bc_outfile_path
