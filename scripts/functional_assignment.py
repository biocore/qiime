#!/usr/bin/env python
# File created on 09 Feb 2010

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso", "Kyle Bittinger",
               "Jens Reeder", "William Walters", "Jose Carlos Clemente Litran"] 
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os.path import splitext, split, exists, abspath
from qiime.util import (make_option, parse_command_line_parameters, create_dir)
from qiime.pick_otus  import BlastxOtuPicker

assignment_constructors = {'blastx':BlastxOtuPicker}

script_info={}
script_info['brief_description'] = """ Script for performing functional assignment of reads against a reference database """
script_info['script_description'] = """ """
script_info['script_usage'] = [("""""","""Run assignment with BLAST using default parameters""","""%prog -i query_nt.fasta -r refseqs_pr.fasta -o blast_assigned_functions/""")]
script_info['script_usage'].append(("""""","""Run assignment with BLAST using scricter e-value threshold""","""%prog -i query_nt.fasta -r refseqs_pr.fasta -o blast_assigned_function_strict/ -e 1e-70"""))
script_info['output_description'] = """ """

script_info['required_options'] = [
    make_option('-i', '--input_seqs_filepath',type='existing_filepath',
        help='Path to input sequences file'),
    ]

script_info['optional_options'] = [
    make_option('-m', '--assignment_method', type='choice',
        choices=assignment_constructors.keys(), default = "blastx",
        help=('Method for picking OTUs.  Valid choices are: ' +\
              ', '.join(assignment_constructors.keys()) +\
              '. [default: %default]')),

    make_option('-o', '--output_dir',type='new_dirpath',\
        help=('Path to store result file '
              '[default: ./<METHOD>_assigned_functions/]')),
              
    make_option('-r', '--refseqs_fp',type='existing_filepath',
        help=('Path to reference sequences to search against when using -m '
              'blast [default: %default]')),
              
    make_option('-b', '--blast_db',type='string',
        help=('Pre-existing database to blast against when using -m blast '
              '[default: %default]')),
              
    make_option('--min_aligned_percent',
        help=('Minimum percent of query sequence that can be aligned to consider a match'
              ' [default: %default]'),default=0.50,type='float'),
              
    make_option('-s', '--min_percent_similarity', type='float', default=0.75,
        help=('Minimum percent similarity to consider a match [default: %default]')),
              
    make_option('-e', '--max_e_value', type='float', default=1e-10,
        help=('Max e-value to consider a match [default: %default]')),
]

script_info['version'] = __version__

def main():
    # Parse the command line parameters
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    
    # Create local copies of the options to avoid repetitive lookups
    assignment_method = opts.assignment_method
    refseqs_fp = opts.refseqs_fp
    if refseqs_fp:
        refseqs_fp = abspath(refseqs_fp)
    blast_db = opts.blast_db
    min_percent_similarity = opts.min_percent_similarity
    min_aligned_percent = opts.min_aligned_percent
    max_e_value = opts.max_e_value
    verbose = opts.verbose
    
    if (assignment_method == 'blastx' and
        refseqs_fp == None and
        blast_db == None):
           option_parser.error('blastx requires refseqs_fp or blast_db')
    
    # use the otu_picking_method value to get the otu picker constructor
    assignment_constructor = assignment_constructors[assignment_method]
    
    # split the input filepath into components used for generating
    # the output file name
    input_seqs_filepath = opts.input_seqs_filepath
    input_seqs_dir, input_seqs_filename = split(input_seqs_filepath)
    input_seqs_basename, ext = splitext(input_seqs_filename)
    
    # create the output directory name (if not provided) and 
    # create it if it doesn't already exist
    output_dir = opts.output_dir or assignment_method + '_assigned_functions'
    create_dir(output_dir, fail_on_exist=False)
    
    # Create the output and log file names
    result_path = '%s/%s_clusters.txt' % (output_dir,input_seqs_basename)
    log_path = '%s/%s_clusters.log' % (output_dir,input_seqs_basename)
    failure_path = '%s/%s_failures.txt' % (output_dir,input_seqs_basename)
    
    if assignment_method == 'blastx':
        params = {'max_e_value':max_e_value,
                  'Similarity': min_percent_similarity,
                  'min_aligned_percent':min_aligned_percent}
        otu_picker = assignment_constructor(params)
        otu_picker(input_seqs_filepath,
                   result_path=result_path, 
                   log_path=log_path,
                   blast_db=opts.blast_db,
                   refseqs_fp=opts.refseqs_fp)
    else:
        ## other -- shouldn't be able to get here as a KeyError would have
        ## been raised earlier
        raise ValueError, "Unknown functional assignment method: %s" % assignment_method

if __name__ == "__main__":
    main()

