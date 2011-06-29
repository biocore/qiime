#!/usr/bin/env python
# File created on 08 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"


from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.util import extract_seqs_by_sample_id
from cogent.parse.fasta import MinimalFastaParser

#extract_seqs_by_sample_id.py
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""Extract sequences based on the SampleID"""
script_info['script_description']="""This script creates a fasta file which will contain only sequences that ARE associated with the supplied sampleIDs(-s), OR all sequences that are NOT associated with the supplied sampleIDs (-n)"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE associated with sample ids S2, S3, \
S4 (-s). As always, sample IDs are case-sensitive:""","""extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4"""))
script_info['script_usage'].append(("""""","""Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences  THAT ARE NOT (-n) associated with sample ids S2, S3, S4 (-s). As always, sample IDs are case-sensitive:""","""extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4 -n"""))
script_info['output_description']="""The script produces a fasta file containing containing only the specified SampleIDs."""
script_info['required_options']=[\
options_lookup['fasta_as_primary_input'],\
make_option('-s','--sample_ids',\
help="comma-separated sample_ids to include in output fasta file"+\
"(or exclude if -n=True)"),\
make_option('-o','--output_fasta_fp',\
help='the output fasta file')]
script_info['optional_options']=[\
make_option('-n','--negate',action='store_true',default=False,\
help='negate the sample ID list (i.e., output sample '+\
'ids not passed via -s) [default: %default]')]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
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
