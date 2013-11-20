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
from cogent.app.fastq_join import run_fastqjoin
from cogent.app.seqprep import run_seqprep
from cogent.app.flash import run_flash
from cogent.app.pandaseq import run_pandaseq

join_method_constructors = {}
join_method_names = {'fastq-join':run_fastqjoin,
                     'SeqPrep':run_seqprep,
                     'flash':run_flash,
                     'pandaseq':run_pandaseq}

def read_bc_to_dict(barcodes):
    """Returns {fastq.header:fastq_data}.
        barcodes : barcodes file handle
    """
    bc_dict = {}
    for data in MinimalFastqParser(barcodes, strict=False):
        curr_label = data[0].strip()
        curr_seq = data[1].strip()
        curr_qual = data[2].strip()
        bc_dict[curr_label] = (curr_seq, curr_qual)
    return bc_dict

def remove_unused_barcodes(assembly, bc_dict, outfile):
    """Writes barcodes to file that only exist in assembly.
       We process this way so that the barcode and assembly fastq 
       files are in the same order!

    assembly : joined paired-ends file-handle
    bc_dict : output from read_bc_to_dict function
    outfile : outfile handle

    """
    for data in MinimalFastqParser(assembly, strict=False):
        curr_label = data[0].strip()
        try:
            curr_seq,curr_qual = bc_dict[curr_label]
            fastq_string = '@%s\n%s\n+\n%s\n' % (curr_label, curr_seq,\
                                                 curr_qual)
            outfile.write(fastq_string)
        except KeyError:
            print "%s not found" % curr_label


