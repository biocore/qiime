#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso, Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from os.path import split, splitext

from qiime.identify_chimeric_seqs import blast_fragments_identify_chimeras

options_lookup = get_options_lookup()

#identify_chimeric_seqs.py
script_info={}
script_info['brief_description']="""Identify chimeric sequences in input FASTA file"""
script_info['script_description']="""A FASTA file of sequences, can be screened to remove chimeras (sequences generated due to the PCR amplification of multiple templates or parent sequences). QIIME currently includes a single taxonomy-assignment-based approach, blast_fragments, for identifying sequences as chimeric. 

The reference sequences (-r) and id-to-taxonomy map (-t) provided are the same format as those provided to assign_taxonomy.py. The reference sequences are in fasta format, and the id-to-taxonomy map contains tab-separated lines where the first field is a sequence identifier, and the second field is the taxonomy separated by semi-colons (e.g., Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium). The reference collection should be derived from a chimera-checked datbase (such as the full greengenes database), and filtered to contain only sequences at, for example, a maximum of 97% sequence identity."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""The blast_fragments chimera detection method is the default method used by identify_chimeric_seqs.py. For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in assign_taxonomy.py) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to the directory "identify_chimeras/" and using default parameters (e.g. chimera detection method ("-m blast_fragments"), number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")):""","""identify_chimeric_seqs.py -i repr_set_seqs.fasta -t taxonomy_assignment.txt -r ref_seq_set.fna -o identify_chimeras/"""))
script_info['output_description']="""The result of identify_chimeric_seqs.py is a text file that identifies which sequences are chimeric."""
script_info['required_options']=[options_lookup['fasta_as_primary_input']]

chimera_detection_method_choices = ['blast_fragments']

script_info['optional_options']=[\
    make_option('-t', '--id_to_taxonomy_fp',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a comma-separated '
         'list. [default: %default; REQUIRED when method is blast_fragments]'),

    make_option('-r', '--reference_seqs_fp',
        help='Path to reference sequences (used to build a blast db).'
        '[default: %default; REQUIRED when method is blast_fragments'+\
         ' if no blast_db is provided]'),
        
    make_option('-b', '--blast_db',
        help='Database to blast against. Must provide either --blast_db or '
        '--reference_seqs_fp when method is blast_fragments [default: %default]'),
        
    make_option('-m','--chimera_detection_method',\
          type='choice',help='Chimera detection method [default:%default]',\
          choices=chimera_detection_method_choices, default='blast_fragments'),
          
    make_option('-n','--num_fragments',\
          type='int',help='Number of fragments to split sequences into' +\
          ' (i.e., number of expected breakpoints + 1) [default: %default]',\
          default=3),
          
    make_option('-d','--taxonomy_depth',\
          type='int',help='Number of taxonomic divisions to consider' +\
          ' when comparing taxonomy assignments [default: %default]',\
          default=4),
          
    make_option('-e','--max_e_value',\
          type='float',help='Max e-value to assign taxonomy' +\
          ' [default: %default]', default=1e-30),
          
    make_option('-o', '--output_fp',
        help='Path to store output [derived from input_seqs_fp]')
]
script_info['version'] = __version__

def main():
    """Run chimera checker with given options>"""

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    #additional option checks
    if opts.chimera_detection_method == 'blast_fragments':
        if not (opts.blast_db or opts.reference_seqs_fp):
            option_parser.error('Must provide either --blast_db or'+\
                ' --reference_seqs_fp and --id_to_taxonomy_fp when'+\
                ' method is blast_fragments.')
        if not opts.id_to_taxonomy_fp:
            option_parser.error('Must provide --id_to_taxonomy_fp when method'+\
                ' is blast_fragments.')

    if opts.num_fragments < 2:
        option_parser.error('Invalid number of fragments (-n %d) Must be >= 2.' \
         % opts.num_fragments)

    verbose = opts.verbose #not used yet ...
    input_seqs_fp = opts.input_fasta_fp
    id_to_taxonomy_fp = opts.id_to_taxonomy_fp
    reference_seqs_fp = opts.reference_seqs_fp
    chimera_detection_method = opts.chimera_detection_method
    num_fragments = opts.num_fragments
    output_fp = opts.output_fp
    taxonomy_depth = opts.taxonomy_depth
    max_e_value = opts.max_e_value
    blast_db = opts.blast_db
    
    if not output_fp:
        input_basename = splitext(split(input_seqs_fp)[1])[0]
        output_fp = '%s_chimeric.txt' % input_basename
    
    if chimera_detection_method == 'blast_fragments':
        blast_fragments_identify_chimeras(input_seqs_fp,
            id_to_taxonomy_fp,\
            reference_seqs_fp,blast_db=blast_db,
            num_fragments=opts.num_fragments,\
            max_e_value=max_e_value,\
            output_fp=output_fp,
            taxonomy_depth=taxonomy_depth)        

if __name__ == "__main__":
    main()
