#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso","Jens Reeder"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from qiime.util import make_option
from cogent.parse.fasta import MinimalFastaParser
from cogent.parse.fastq import MinimalFastqParser
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import fields_to_dict
from qiime.filter import (filter_fasta, filter_fastq,
                          get_seqs_to_keep_lookup_from_seq_id_file,
                          get_seqs_to_keep_lookup_from_fasta_file,
                          sample_ids_from_metadata_description)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "This script can be applied to remove sequences from a fasta or fastq file based on input criteria."
script_info['script_description'] = ""
script_info['script_usage'] = [
 ("Keep all sequences that show up in an OTU map","",
 "filter_fasta.py -f inseqs.fasta -o filtered_seqs.fasta -m uclust_ref_otus.txt"),
 ("Discard all sequences that show up in chimera checking output. NOTE: It is very important to pass -n here as this tells the script to negate the request, or discard all sequences that are listed via -s. This is necessary to remove the identified chimeras from inseqs.fasta","",
 "filter_fasta.py -f inseqs.fasta -o non_chimeric_seqs.fasta -s chimeric_seqs.txt -n"),
 ("Keep all sequences from as fasta file that are listed in a text file","",
 "filter_fasta.py -f inseqs.fasta -o filtered_seqs.fasta -s seqs_to_keep.txt"),
 ("Keep all sequences from a fastq file that are listed in a text file (note: file name must end with .fastq to support fastq filtering)","",
 "filter_fasta.py -f inseqs.fastq -o filtered_seqs.fasta -s seqs_to_keep.txt")]
script_info['output_description']= ""
script_info['required_options'] = [\
 options_lookup['input_fasta'],
 make_option('-o','--output_fasta_fp',help='the output fasta filepath')
]
script_info['optional_options'] = [\
 make_option('-m','--otu_map',
  help='an OTU map where sequences ids are those which should be retained'),\
 make_option('-s','--seq_id_fp', 
  help='A list of sequence identifiers (or tab-delimited lines with'
  ' a seq identifier in the first field) which should be retained'),\
 make_option('-a','--subject_fasta_fp',
  help='A fasta file where the seq ids should be retained.'),\
 make_option('-p','--seq_id_prefix',
  help='keep seqs where seq_id starts with this prefix'),\
 make_option('-n','--negate', help='discard passed seq ids rather than'
  ' keep passed seq ids [default: %default]', default=False, 
  action='store_true'),
 make_option('--mapping_fp',
  help='mapping file path (for use with --valid_states) [default: %default]'),
 make_option('--valid_states',
  help='description of sample ids to retain (for use with --mapping_fp) [default: %default]')
]
script_info['version'] = __version__

def filter_fasta_fp(input_seqs_fp,output_seqs_fp,seqs_to_keep,negate=False):
    """Filter a fasta file to include only sequences listed in seqs_to_keep """
    input_seqs = MinimalFastqParser(open(input_seqs_fp,'U'))
    output_f = open(output_seqs_fp,'w')
    return filter_fasta(input_seqs,output_f,seqs_to_keep,negate)

def filter_fastq_fp(input_seqs_fp,output_seqs_fp,seqs_to_keep,negate=False):
    """Filter a fastq file to include only sequences listed in seqs_to_keep """
    input_seqs = MinimalFastqParser(open(input_seqs_fp,'U'),strict=False)
    output_f = open(output_seqs_fp,'w')
    return filter_fastq(input_seqs,output_f,seqs_to_keep,negate)

def get_seqs_to_keep_lookup_from_otu_map(seqs_to_keep_f):
    """Generate a lookup dictionary from an OTU map"""
    otu_map = fields_to_dict(seqs_to_keep_f)
    seqs_to_keep = []
    for seq_ids in otu_map.values():
        seqs_to_keep += seq_ids
    return {}.fromkeys(seqs_to_keep)

def get_seqs_to_keep_lookup_from_prefix(fasta_f,prefix):
    seqs_to_keep = [seq_id
                    for seq_id, seq in MinimalFastaParser(fasta_f)
                    if seq_id.startswith(prefix)]
    return {}.fromkeys(seqs_to_keep)

def get_seqs_to_keep_lookup_from_mapping_file(fasta_f,mapping_f,valid_states):
    sample_ids = {}.fromkeys(\
     sample_ids_from_metadata_description(mapping_f,valid_states))
    seqs_to_keep = []
    for seq_id, seq in MinimalFastaParser(fasta_f):
        if seq_id.split('_')[0] in sample_ids:
            seqs_to_keep.append(seq_id)
        else:
            continue
    return {}.fromkeys(seqs_to_keep)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    negate = opts.negate
    error_msg = "Must pass exactly one of -a, -s, -p, -m, or --valid_states and --mapping_fp."
    if 1 != sum(map(bool,[opts.otu_map,
                          opts.seq_id_fp,
                          opts.subject_fasta_fp,
                          opts.seq_id_prefix,
                          opts.mapping_fp and opts.valid_states])): 
        option_parser.error(error_msg)

    if opts.otu_map:
        seqs_to_keep_lookup =\
         get_seqs_to_keep_lookup_from_otu_map(
         open(opts.otu_map,'U'))
    elif opts.seq_id_fp:
        seqs_to_keep_lookup =\
         get_seqs_to_keep_lookup_from_seq_id_file(
         open(opts.seq_id_fp,'U'))
    elif opts.subject_fasta_fp:
        seqs_to_keep_lookup =\
         get_seqs_to_keep_lookup_from_fasta_file(
         open(opts.subject_fasta_fp,'U'))
    elif opts.seq_id_prefix:
        seqs_to_keep_lookup =\
         get_seqs_to_keep_lookup_from_prefix(
         open(opts.input_fasta_fp),opts.seq_id_prefix)
    elif opts.mapping_fp and opts.valid_states:
        seqs_to_keep_lookup =\
         get_seqs_to_keep_lookup_from_mapping_file(
          open(opts.input_fasta_fp,'U'),
          open(opts.mapping_fp,'U'),
          opts.valid_states)
    else:
        option_parser.error(error_msg)
    
    if opts.input_fasta_fp.endswith('.fastq'):
        filter_fp_f = filter_fastq_fp
    else:
        filter_fp_f = filter_fasta_fp
    
    filter_fp_f(opts.input_fasta_fp,
                opts.output_fasta_fp,
                seqs_to_keep_lookup,
                negate)

if __name__ == "__main__":
    main()
