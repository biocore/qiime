#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso, Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald","Jens Reeder",
                "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from os.path import split, splitext

from qiime.identify_chimeric_seqs import blast_fragments_identify_chimeras,\
    chimeraSlayer_identify_chimeras

options_lookup = get_options_lookup()

#identify_chimeric_seqs.py
script_info={}
script_info['brief_description']="""Identify chimeric sequences in input FASTA file"""
script_info['script_description']="""A FASTA file of sequences, can be screened to remove chimeras (sequences generated due to the PCR amplification of multiple templates or parent sequences). QIIME currently includes a taxonomy-assignment-based approach, blast_fragments, for identifying sequences as chimeric and the ChimeraSlayer algorithm. 

1. Blast_fragments approach: 

The reference sequences (-r) and id-to-taxonomy map (-t) provided are the same format as those provided to assign_taxonomy.py. The reference sequences are in fasta format, and the id-to-taxonomy map contains tab-separated lines where the first field is a sequence identifier, and the second field is the taxonomy separated by semi-colons (e.g., Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium). The reference collection should be derived from a chimera-checked database (such as the full greengenes database), and filtered to contain only sequences at, for example, a maximum of 97% sequence identity.

2. ChimeraSlayer:

ChimeraSlayer uses BLAST to identify potential chimera parents and computes the optimal branching alignment of the query against two parents.
We suggest to use the pynast aligned representative sequences as input.
"""

script_info['script_usage']=[]
script_info['script_usage'].append(("""blast_fragments example""","""For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in assign_taxonomy.py) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to the directory "identify_chimeras/" and using default parameters (e.g. chimera detection method ("-m blast_fragments"), number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")):""","""%prog -i repr_set_seqs.fasta -t taxonomy_assignment.txt -r ref_seq_set.fna -m blast_fragments -o chimeric_seqs_blast.txt"""))

script_info['script_usage'].append(("""ChimeraSlayer Example:""","""Identify chimeric sequences using the ChimeraSlayer algorithm against a user provided reference data base. The input sequences need to be provided in aligned (Py)Nast format. The reference data base needs to be provided as aligned FASTA (-a). Note that the reference database needs to be the same that was used to build the alignment of the input sequences!""",
                                    """%prog -m ChimeraSlayer -i repr_set_seqs_aligned.fasta -a ref_seq_set_aligned.fasta -o chimeric_seqs_cs.txt"""))

script_info['output_description']="""The result of identify_chimeric_seqs.py is a text file that identifies which sequences are chimeric."""
script_info['required_options']=[options_lookup['fasta_as_primary_input']]

chimera_detection_method_choices = ['blast_fragments','ChimeraSlayer']

script_info['optional_options']=[\
    make_option('-t', '--id_to_taxonomy_fp',type='existing_filepath',
        help='Path to tab-delimited file mapping sequences to assigned '
         'taxonomy. Each assigned taxonomy is provided as a comma-separated '
         'list. [default: %default; REQUIRED when method is blast_fragments]'),

    make_option('-r', '--reference_seqs_fp',type='existing_filepath',
        help='Path to reference sequences (used to build a blast db when method blast_fragments). '
        '[default: %default; REQUIRED when method blast_fragments'+\
         ' if no blast_db is provided;]'),

    make_option('-a', '--aligned_reference_seqs_fp',type='existing_filepath',
        help='Path to (Py)Nast aligned reference sequences. '
        'REQUIRED when method ChimeraSlayer [default: %default]'),       

    make_option('-b', '--blast_db',type='blast_db',
        help='Database to blast against. Must provide either --blast_db or '
        '--reference_seqs_fp when method is blast_fragments [default: %default]'),
        
    make_option('-m','--chimera_detection_method',\
          type='choice',help='Chimera detection method. Choices: '+\
                    " or ".join(chimera_detection_method_choices) +\
                    '. [default:%default]',\
          choices=chimera_detection_method_choices, default='ChimeraSlayer'),
          
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

    make_option('-R','--min_div_ratio',\
                    type='float',help='min divergence ratio '+\
                    '(passed to ChimeraSlayer). If set to None uses ' +\
                    'ChimeraSlayer default value. '+\
                    ' [default: %default]', default=None),       
  
    make_option('-k','--keep_intermediates',\
                    action='store_true',help='Keep intermediate files, '+\
                    'useful for debugging ' +\
                    ' [default: %default]', default=False),  
    
    make_option('-o', '--output_fp',type='new_filepath',
        help='Path to store output [default: derived from input_seqs_fp]')
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

    elif opts.chimera_detection_method == 'ChimeraSlayer':
        if not opts.aligned_reference_seqs_fp:
            option_parser.error("Must provide --aligned_reference_seqs_fp "+\
                                    "when using method ChimeraSlayer")
            
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
    keep_intermediates = opts.keep_intermediates
    
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

    elif chimera_detection_method == 'ChimeraSlayer':
         chimeraSlayer_identify_chimeras(input_seqs_fp, 
                                         output_fp=output_fp,
                                         db_FASTA_fp=opts.reference_seqs_fp,
                                         db_NAST_fp=opts.aligned_reference_seqs_fp,
                                         min_div_ratio=opts.min_div_ratio,
                                         keep_intermediates=keep_intermediates)

if __name__ == "__main__":
    main()
