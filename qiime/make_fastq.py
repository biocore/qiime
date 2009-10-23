#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"


"""Generates fastq file for SRA submission from paired fasta and qual files.

The fastq format for each record is as follows:
@seq_id [and optional description]
seq as bases
+[optionally with repeat of seq_id and repeat line]
qual scores as string of chr(33+qual)

The SRA currently requires a separate fastq file for each library, split by
library id. This code takes the output from split_libraries.py and the
corresponding qual files, pulls the qual info by id, and writes everything
either to one file or to per-library files.
"""
from qiime.split_libraries import qual_scores
from cogent.parse.fasta import MinimalFastaParser
from os import mkdir
from collections import defaultdict
from optparse import OptionParser

def make_fastq_rec(header, seq, qual, offset=33):
    """Makes fasta string given header, seq, qual scores as integers"""
    result = []
    if header.startswith('>'):
        header = header[1:]
    result.append('@'+header)
    result.append(seq)
    result.append('+'+header)
    result.append(''.join(map(chr, [33+i for i in qual])))
    return '\n'.join(result)

def split_lib_transform(header):
    """Transform header of split_libraries.py for submission.

    Specifically:
    renames orig_bc field to 'barcode'
    returns machine id to retrieve qual file

    Result is tuple of (new_header, qual_id)
    """
    if header.startswith('>'):
        header = header[1:]
    fields = header.split()
    lib_id = fields[0]
    qual_id = fields[1]
    bc = fields[2].split('=')[1]
    return ' '.join([lib_id, 'read_id='+qual_id, 'barcode='+bc]), qual_id

def iter_fastq(in_fasta, quals, label_transform=split_lib_transform):
    """Iterate over fastq records, yields seq id of each"""
    for label, seq in MinimalFastaParser(in_fasta):
        new_label, qual_id = label_transform(label)
        seq_id = label.split()[0]
        if seq_id.startswith('>'):
            seq_id = seq_id[1:]
        qual = quals[qual_id]
        yield make_fastq_rec(new_label, seq, qual), seq_id

def make_fastq_single(in_fasta, quals, out_fp, \
    label_transform=split_lib_transform):
    """Makes a single fastq file with all the data"""
    outfile = open(out_fp, 'w')
    for rec, seq_id in iter_fastq(in_fasta, quals, label_transform):
        outfile.write(rec+'\n')
    outfile.close()

def make_fastq_multi(in_fasta, quals, out_fp, \
    label_transform=split_lib_transform):
    """Makes a fastq file for each library (assumed to be libid_seqid in label)
    
    WARNING: this implementation reads everything into memory, so won't scale.
    """
    mkdir(out_fp)
    seen_libs = defaultdict(list)
    for rec, label in iter_fastq(in_fasta, quals, label_transform):
        lib_id, seq_id = label.rsplit('_',1)
        seen_libs[lib_id].append(rec)
    for lib, recs in seen_libs.items():
        if lib is None: #skip the seqs we couldn't assign to a library
            continue
        outfile = open(out_fp+'/'+lib+'.fastq', 'w')
        outfile.write('\n'.join(recs))
        outfile.close()

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage =\
     'usage: %prog [options] -f input_fasta_filepath -q input_qual_filepath'
    version = 'Version: %prog ' +  __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-f', '--fasta', dest='fasta_fp',
        help='name of fasta output of split_libraries.py [Required]')

    parser.add_option('-q', '--qual', dest='qual_fps',
        help='names of qual files, comma-delimited [Required]')

    parser.add_option('-o','--result_fp',action='store',\
          type='string',dest='result_fp',help='Path to store '+\
          'results [default: <input_sequences_filename>.fastq]')

    parser.add_option('-s', '--split', dest='split', action='store_true',
        help='make separate file for each library [default:False]', 
        default=False)
          
    parser.set_defaults(result_fp=None)

    opts,args = parser.parse_args()
    return opts,args

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    #verbose = opts.verbose
    in_fasta = open(opts.fasta_fp, 'U')
    quals = qual_scores([open(f, 'U') for f in opts.qual_fps.split(',')])
    if not opts.result_fp:
        opts.result_fp = opts.fasta_fp + '.fastq'

    if opts.split:
        make_fastq_multi(in_fasta, quals, opts.result_fp)
    else:
        make_fastq_single(in_fasta, quals, opts.result_fp)


