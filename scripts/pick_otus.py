#!/usr/bin/env python
# File created on 09 Feb 2010

__author__ = "William Walters"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger","Jens Reeder","William Walters"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.pick_otus  import otu_picking_method_constructors,\
 otu_picking_method_choices, MothurOtuPicker
from os.path import splitext, split
from os import makedirs

script_description = """ This module creates clusters of OTUs from a given 
input fasta file, which are written to an output .txt file name based on the
fasta file name.  By default, the OTU picking method is cdhit, and the
percent similarity for clusters is 0.97. """

script_usage = """
    # Get detailed usage information
    python pick_otus.py -h
    
    # Pick OTUs from at_inseqs.fasta (-i) with cdhit (default) and
    # store output files in ./cdhit_picked_otus/ (default).
    python pick_otus.py -i at_inseqs.fasta

    # Pick OTUs from at_inseqs.fasta with Mothur, using a nearest-
    # neighbor clustering algorithm and a similarity threshold of 90%.
    # Output files will be stored in ./mothur_picked_otus/ (default).
    python pick_otus.py -i at_inseqs.fasta -m mothur -c nearest -s 0.90
    
    # Pick OTUs from inseqs.fasta (-i) using the BLAST otu picker (-m)
    # using ref_seqs.fasta to build a blast database on-the-fly. (Note that
    # a pre-existing blast database can also be provided via the -b parameter).
    python pick_otus.py -m blast -i inseqs.fasta -r ref_seqs.fasta"""
    
required_options = [\
 make_option('-i','--input_seqs_filepath', help='Path to input sequences file')
]


optional_options = [\
 make_option('-m', '--otu_picking_method', type='choice',
        choices=otu_picking_method_choices, default = "cdhit",
        help=('Method for picking OTUs.  Valid choices are: ' +\
        ', '.join(otu_picking_method_choices) +\
        '. The mothur method requires an input file ' +\
        'of aligned sequences [default: %default]')),\
 make_option('-c', '--clustering_algorithm', type='choice',
        choices=MothurOtuPicker.ClusteringAlgorithms, default = "furthest",
        help=('Clustering algorithm for mothur otu picking method.  Valid '
        'choices are: nearest, furthest, and average [default: %default]')),\
 make_option('-M','--max_cdhit_memory',type=int,\
          help='max available memory to cdhit (cd-hit\'s -M)'+\
          ' (Mbyte) [default: %default]',default=400),\
 make_option('-o','--output_dir',\
          help='Path to store '+\
          'result file [default: ./<OTU_METHOD>_picked_otus/]'),\
 make_option('-r','--refseqs_fp',
          help='Path to reference sequences to blast against when'+\
          ' using -m blast [default: %default]'),\
 make_option('-b','--blast_db',
          help='Pre-existing database to blast against when'+\
          ' using -m blast [default: %default]'),\
 make_option('-s','--similarity',action='store', default = 0.97,\
          type='float',dest='similarity',help='Sequence similarity '+\
          'threshold (for cdhit or uclust) [default: %default]'),\
 make_option('-e','--max_e_value',action='store', default = 1e-10,\
          type='float',dest='max_e_value',help='Max E-value when '+\
          'clustering with BLAST [default: %default]'),\
 make_option('-q','--trie_reverse_seqs',action='store_true', default = False,\
          help='Reverse seqs before picking OTUs with the Trie OTU'+\
          ' picker for suffix (rather than prefix) collapsing'+\
          ' [default: %default]'),\
 make_option('-n','--prefix_prefilter_length',\
          type=int,help='prefilter data so seqs with identical first '+\
          'prefix_prefilter_length are automatically grouped into a '+\
          'single OTU; useful for large sequence collections where OTU '+\
          'picking doesn\'t scale well '+\
          '[default: %default; 100 is a good value]', default = 50),\
 make_option('-t','--trie_prefilter',\
          action='store_true', default = False,\
          help='prefilter data so seqs which are identical prefixes '+\
          'of a longer seq are automatically grouped into a '+\
          'single OTU; useful for large sequence collections where OTU '+\
          'picking doesn\'t scale well '+\
          '[default: %default]'),\
 make_option('-p','--prefix_length', default = 50,\
          type=int,help='prefix length when using the prefix_suffix'+\
          ' otu picker; WARNING: CURRENTLY DIFFERENT FROM'+\
          ' prefix_prefilter_length (-n)! [default: %default]'),\
 make_option('-u','--suffix_length', default = 50,\
          type=int,help='suffix length when using the prefix_suffix'+\
          ' otu picker [default: %default]') ]
          
def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    prefix_prefilter_length = opts.prefix_prefilter_length
    otu_picking_method = opts.otu_picking_method
    prefix_length = opts.prefix_length
    suffix_length = opts.suffix_length
    trie_prefilter = opts.trie_prefilter
    trie_reverse_seqs = opts.trie_reverse_seqs
    refseqs_fp = opts.refseqs_fp
    blast_db = opts.blast_db
    similarity = opts.similarity
    
    if otu_picking_method == 'cdhit' and similarity < 0.80:
        option_parser.error('cdhit requires similarity >= 0.80.')
            
    if otu_picking_method == 'blast' and \
       refseqs_fp == None and \
       blast_db == None:
           option_parser.error('blast requires refseqs_fp or blast_db')
 
    otu_picker_constructor =\
     otu_picking_method_constructors[otu_picking_method]
     
    input_seqs_filepath = opts.input_seqs_filepath
    input_seqs_dir, input_seqs_filename = split(input_seqs_filepath)
    input_seqs_basename, ext = splitext(input_seqs_filename)
    
    output_dir = opts.output_dir or otu_picking_method + '_picked_otus'
    try:
        makedirs(output_dir)
    except OSError:
        # output_dir exists
        pass
        
    result_path = '%s/%s_otus.txt' % (output_dir,input_seqs_basename)
    log_path = '%s/%s_otus.log' % (output_dir,input_seqs_basename)
    
    if otu_picking_method == 'cdhit':
        params = {'Similarity':opts.similarity,\
                  '-M':opts.max_cdhit_memory}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,\
         result_path=result_path,log_path=log_path,\
         prefix_prefilter_length=prefix_prefilter_length,\
         trie_prefilter=trie_prefilter)
    elif otu_picking_method == 'uclust':
        params = {'Similarity':opts.similarity}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,\
         result_path=result_path,log_path=log_path,\
         prefix_prefilter_length=prefix_prefilter_length,\
         trie_prefilter=trie_prefilter)
    elif otu_picking_method == 'prefix_suffix':
        otu_picker = otu_picker_constructor({})
        otu_picker(input_seqs_filepath,\
         result_path=result_path,log_path=log_path,\
         prefix_length=prefix_length,suffix_length=suffix_length)
    elif otu_picking_method == 'mothur':
        params = {'Similarity': opts.similarity,
                  'Algorithm': opts.clustering_algorithm}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)
    elif otu_picking_method == 'trie':
        params = {'Reverse':trie_reverse_seqs}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path)
    elif otu_picking_method == 'blast':
        params = {'max_e_value':opts.max_e_value,\
                  'Similarity': opts.similarity}
        otu_picker = otu_picker_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, log_path=log_path,\
                   blast_db=opts.blast_db,refseqs_fp=opts.refseqs_fp)


if __name__ == "__main__":
    main()

