#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger", "Antonio Gonzalez Pena"] 
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from os import system, remove, path, mkdir
from os.path import split, splitext
from qiime.assign_taxonomy import BlastTaxonAssigner, RdpTaxonAssigner

script_description = """Contains code for assigning taxonomy, using several techniques.

This module has the responsibility for taking a set of sequences and
providing a taxon assignment for each sequence."""

script_usage = """ Assign taxonomy of sequences in inseqs.fasta (-i) using the RDP classifier 
 (-m). Output files will be written to rdp_assigned_taxonomy (default).
 assign_taxonomy.py -i inseqs.fasta -m rdp
 
 Assign taxonomy of sequences in inseqs.fasta (-i) using the RDP classifier
 (-m) trained on-the-fly from provided refseqs and taxon assignments (-r, -t) 
 respectively. Output files will be written to custom_rdp.
 assign_taxonomy.py -i inseqs.fasta -r refseqs.fasta -t id_to_taxonomy.txt -m rdp -o custom_rdp
 
 Assign taxonomy of sequences in inseqs.fasta (-i) using BLAST
 (-m) against provided refseqs and taxon assignments (-r, -t) 
 respectively. Output files will be written to blast_assigned_taxonomy (default).
 assign_taxonomy.py -i inseqs.fasta -r at_refseqs.fasta -t at_id_to_taxonomy.txt -m blast
"""

assignment_method_constructors = {
    'blast': BlastTaxonAssigner,
    'rdp': RdpTaxonAssigner,
    }
 
assignment_method_choices = assignment_method_constructors.keys()
    
required_options = [\
 make_option('-i', '--input_seqs_fp',
        help='Path to fasta file of sequences to be assigned [REQUIRED]')
]

optional_options = [\
 make_option('-t', '--id_to_taxonomy_fp',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a comma-separated '
         'list. For assignment with rdp, each assigned taxonomy must be '
         'exactly 6 levels deep. [default: %default; REQUIRED when method is '
         'blast]'),\
 make_option('-r', '--reference_seqs_fp',
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp, they '
        'are used as training sequences for the classifier.'
        '[default: %default; REQUIRED if -b is not provided when method is blast]'),\
 make_option('-m','--assignment_method',\
          type='choice',help='Taxon assignment method [default:%default]',\
          choices=assignment_method_choices, default="rdp"),\
 make_option('-b', '--blast_db',
        help='Database to blast against.  Must provide either --blast_db or '
        '--reference_seqs_db for assignment with blast [default: %default]'),\
 make_option('-c', '--confidence', type='float',
        help='Minimum confidence to record an assignment, only used for rdp '
        'method [default: %default]', default=0.80),\
 make_option('-e', '--e_value', type='float',
        help='Maximum e-value to record an assignment, only used for blast '
        'method [default: %default]',default=0.001),\
 make_option('-o','--output_dir',\
          help='Path to store result file '+\
          '[default: <ASSIGNMENT_METHOD>_assigned_taxonomy]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    
    if opts.assignment_method == 'blast':
        if not opts.id_to_taxonomy_fp:
            option_parser.error('Option --id_to_taxonomy_fp is required when ' 
                         'assigning with blast.')
        if not (opts.reference_seqs_fp or opts.blast_db):
            option_parser.error('Either a blast db (via -b) or a collection of '
                         'reference sequences (via -r) must be passed to '
                         'assign taxonomy using blast.')

    if opts.assignment_method == 'rdp':
        if opts.id_to_taxonomy_fp:
            if opts.reference_seqs_fp is None:
                option_parser.error('A filepath for reference sequences must be '
                             'specified (via -r) along with the id_to_taxonomy '
                             'file to train the Rdp Classifier.')
        elif opts.reference_seqs_fp:
                option_parser.error('A filepath for an id to taxonomy map must be '
                             'specified (via -t) along with the reference '
                             'sequences fp to train the Rdp Classifier.')
        else:
            pass

            
    assignment_method = opts.assignment_method
    taxon_assigner_constructor =\
     assignment_method_constructors[assignment_method]
    input_sequences_filepath = opts.input_seqs_fp
    
    try:
        id_to_taxonomy_fp = opts.id_to_taxonomy_fp
        params = {'id_to_taxonomy_filepath':id_to_taxonomy_fp}
    except IndexError:
        params = {}
    
    # Build the output filenames
    output_dir = opts.output_dir or assignment_method + '_assigned_taxonomy'
    try:
        mkdir(output_dir)
    except OSError:
        # output_dir already exists
        pass
        
    fpath, ext = splitext(input_sequences_filepath)
    input_dir, fname = split(fpath)
    result_path = output_dir + '/' + fname + '_tax_assignments.txt'
    log_path = output_dir + '/' + fname + '_tax_assignments.log'
    
    if opts.assignment_method == 'blast':
        # one of these must have a value, otherwise we'd have 
        # an optparse error
        if opts.blast_db:
            params['blast_db'] = opts.blast_db
        else:
            params['reference_seqs_filepath'] = opts.reference_seqs_fp
        params['Max E value'] = opts.e_value

    elif opts.assignment_method == 'rdp':
        params['Confidence'] = opts.confidence
        params['id_to_taxonomy_fp'] = opts.id_to_taxonomy_fp
        params['reference_sequences_fp'] = opts.reference_seqs_fp

    else:
        # should not be able to get here as an unknown classifier would
        # have raised an optparse error
        exit(1)

    taxon_assigner = taxon_assigner_constructor(params)
    taxon_assigner(input_sequences_filepath,\
     result_path=result_path,log_path=log_path)


if __name__ == "__main__":
    main()