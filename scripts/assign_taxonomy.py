#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Greg Caporaso", "Kyle Bittinger", "Antonio Gonzalez Pena"] 
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from os import system, remove, path, mkdir
from os.path import split, splitext
from qiime.assign_taxonomy import BlastTaxonAssigner, RdpTaxonAssigner

assignment_method_constructors = {
    'blast': BlastTaxonAssigner,
    'rdp': RdpTaxonAssigner,
    }
 
assignment_method_choices = assignment_method_constructors.keys()

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Assign taxonomy to each sequence"""
script_info['script_description']="""Contains code for assigning taxonomy, using several techniques.

Given a set of sequences, assign_taxonomy attempts to assign the taxonomy of each sequence. Currently there are two methods implemented: assignment with BLAST and assignment with the RDP classifier. The output of this step is a mapping of input sequence identifiers (1st column of output file) to taxonomy (2nd column) and quality score (3rd column). The sequence identifier of the best BLAST hit is also included if the blast method is used (4th column). """
script_info['script_usage']=[]
script_info['script_usage'].append(("""""","""
Example of consensus lineage: 

The OTU containing 5 sequences annotated as shown below would be assigned to the "Desulfovibrionaceae" level because only 80% of sequences agree with the "LE30" annotation.

* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; 
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30

Assignments are provided in a two column tab-delimited format, which maps input sequence identifiers to assignments. Each assignment is specified as a list of taxa separated by a ';' character.

Example of an assignment output file:

======== =================================================================
AY800210 Archaea;Euryarchaeota;Halobacteriales;uncultured 
EU883771 Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.
EF503699 Archaea;Crenarchaeota;uncultured;uncultured 
DQ260310 Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium 
EF503697 Archaea;Crenarchaeota;uncultured;uncultured
======== =================================================================
""",""""""))
script_info['script_usage'].append(("""Sample Assignment with BLAST:""","""
Taxonomy assignments are made by searching input sequences against a blast database of pre-assigned reference sequences. If a satisfactory match is found, the reference assignment is given to the input sequence. This method does not take the hierarchical structure of the taxonomy into account, but it is very fast and flexible. If a file of reference sequences is provided, a temporary blast database is built on-the-fly. The quality scores assigned by the BLAST taxonomy assigner are e-values.

To assign the sequences to the representative sequence set, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "blast_assigned_taxonomy", you can run the following command:""","""assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt"""))
script_info['script_usage'].append(("""""","""Optionally, the user could changed the E-value ("-e"), using the following command:""","""assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -e 0.01"""))
script_info['script_usage'].append(("""Assignment with the RDP Classifier:""","""The RDP Classifier program (Wang, Garrity, Tiedje, & Cole, 2007) assigns taxonomies by matching sequence segments of length 8 to a database of previously assigned sequences. It uses a naive bayesian algorithm, which means that for each potential assignment, it attempts to calculate the probability of the observed matches, assuming that the assignment is correct and that the sequence segments are completely independent. The RDP Classifier is distributed with a pre-built database of assigned sequence, which is used by default. The quality scores provided by the RDP classifier are confidence values.

To assign the representative sequence set, where the output directory is "rdp_assigned_taxonomy", the you can run the following command:
""","""assign_taxonomy.py -i repr_set_seqs.fasta -m rdp"""))
script_info['script_usage'].append(("""""","""Alternatively, the user could change the minimum confidence score ("-c"), using the following command:""","""assign_taxonomy.py -i repr_set_seqs.fasta -m rdp -c 0.85"""))
script_info['script_usage'].append(("""""","""Note: If a reference set of sequences and taxonomy to id assignment file are provided, the script will use them to generate a new training dataset for the RDP Classifier on-the-fly. Due to limitations in the generation of a training set, each provided assignment must contain exactly 6 taxa in the following order: domain (level=2), phylum (level=3), class (level=4), order (5), family (level=6), and genus (level=7). Additionally, each genus name must be unique, due to the internal algorithm used by the RDP Classifier.
""",""""""))
script_info['output_description']="""The consensus taxonomy assignment implemented here is the most detailed lineage description shared by 90% or more of the sequences within the OTU (this level of agreement can be adjusted by the user). The full lineage information for each sequence is one of the output files of the analysis. In addition, a conflict file records cases in which a phylum-level taxonomy assignment disagreement exists within an OTU (such instances are rare and can reflect sequence misclassification within the greengenes database)."""
script_info['required_options']=[\
   options_lookup['fasta_as_primary_input']\
]
script_info['optional_options']=[\
 make_option('-t', '--id_to_taxonomy_fp',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a semicolon-separated'
         ' list. For assignment with rdp, each assigned taxonomy must be '
         'exactly 6 levels deep. [default: %default; REQUIRED when method is '
         'blast]'),\
 make_option('-r', '--reference_seqs_fp',
        help='Path to reference sequences.  For assignment with blast, these '
        'are used to generate a blast database. For assignment with rdp, they '
        'are used as training sequences for the classifier.'
        '[default: %default; REQUIRED if -b is not provided when method is blast]'),\
 make_option('-p', '--training_data_properties_fp',
        help='Path to ".properties" file in pre-compiled training data for the '
        'RDP Classifier.  This option is overridden by the -t and -r options. '
        '[default: %default]'),\
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
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
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
    input_sequences_filepath = opts.input_fasta_fp
    
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
        params['training_data_properties_fp'] = opts.training_data_properties_fp
    else:
        # should not be able to get here as an unknown classifier would
        # have raised an optparse error
        exit(1)

    taxon_assigner = taxon_assigner_constructor(params)
    taxon_assigner(input_sequences_filepath,\
     result_path=result_path,log_path=log_path)


if __name__ == "__main__":
    main()
