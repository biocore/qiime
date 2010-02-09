#!/usr/bin/env python
# File created on 08 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"


from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.util import extract_seqs_by_sample_id
from cogent.parse.fasta import MinimalFastaParser

script_description = """This script creates a fasta file which will contain \
only sequences that ARE associated with the supplied sampleIDs(-s), OR all sequences \
that are NOT associated with the supplied sampleIDs(-n)"""

script_usage = """Usage: %prog [options] {-i INPUT_FASTA_FP -o OUTPUT_FASTA_FP -s SAMPLE_IDS}

Example usage:

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta \
(-i) containing only the sequences THAT ARE associated with sample ids S2, S3, \
S4 (-s). (As always, sample IDs are case-sensitive.)

python extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta \
(-i) containing only the sequences  THAT ARE NOT (-n) associated with sample \
ids S2, S3, S4 (-s). (As always, sample IDs are case-sensitive.)

python extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4 -n

"""

required_options = [\
make_option('-i','--input_fasta_fp',\
help='the input fasta file'),\
make_option('-s','--sample_ids',\
help="comma-separated sample_ids to include in output fasta file"+\
"(or exclude if -n=True)"),\
make_option('-o','--output_fasta_fp',\
help='the output fasta file')
]

optional_options = [\
make_option('-n','--negate',action='store_true',default='false',\
help='negate the sample ID list (i.e., output sample '+\
'ids not passed via -s) [default: %default]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    
    sample_ids = opts.sample_ids.split(',')
    negate = opts.negate
    
    try:
        seqs = MinimalFastaParser(open(opts.input_fasta_fp))
    except IOError:
        option_parser.error('Cannot open %s. Does it exist? Do you have read access?'%\
         opts.input_fasta_fp)
        exit(1)
        
    try:
        output_fasta_f = open(opts.output_fasta_fp,'w')
    except IOError:
        option_parser.error("Cannot open %s. Does path exist? Do you have write access?" %\
         opts.output_fasta_fp)
        exit(1)
    
    for r in extract_seqs_by_sample_id(seqs,sample_ids,negate):
        output_fasta_f.write('>%s\n%s\n' % r)
    output_fasta_f.close()

if __name__ == "__main__":
    main()