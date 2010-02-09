#!/usr/bin/env python
# File created on 20 Dec 2009.
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
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import qiime_blast_seqs

script_description = """This script is a functionality-limited interface to the \
qiime.util.qiime_blast_seqs function, primarily useful for testing purposes. \
Once that function has been integrated into qiime as the primary blast interface \
it will move to PyCogent. An expanded version of this command line interface may \
replace the script functionality of cogent.app.blast at that point."""

script_usage = """Usage: %prog [options] {-i INPUT_FASTA_FP -r REFSEQS_FP}

Example usage:
Blast all sequences in inseqs.fasta (-i) against a BLAST db constructed \
from refseqs.fasta (-r).

python blast_wrapper.py -i inseqs.fasta -r refseqs.fasta
"""

required_options = [\
make_option('-i','--input_fasta_fp',\
    help='paths to sequences to blast as a fasta file'),\
make_option('-r','--refseqs_fp',\
    help='path to blast database as a fasta file')
]

optional_options = [\
make_option('-n','--num_seqs_per_blast_run', type='int', default='1000', \
help = 'number of sequences passed to each blast call '+\
"- useful for very large sequence collections [default: %default]")
]

def main():
    option_parser, options, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    
    blast_results = qiime_blast_seqs(\
     seqs=MinimalFastaParser(open(options.input_fasta_fp)),\
     refseqs_fp=options.refseqs_fp,\
     seqs_per_blast_run=options.num_seqs_per_blast_run)
     
    for query_id, blast_result in blast_results.items():
        first_blast_result = blast_result[0][0]
        print '%s: %s %s %s' % (\
         query_id,
         first_blast_result['SUBJECT ID'],
         first_blast_result['E-VALUE'],
         first_blast_result['% IDENTITY'])

if __name__ == "__main__":
    main()
    