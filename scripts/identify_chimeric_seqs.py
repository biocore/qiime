#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso, Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from os.path import split, splitext

from qiime.identify_chimeric_seqs import blast_fragments_identify_chimeras

script_description = """Identify chimeric sequences in input FASTA file."""

script_usage = """
 Test whether each seq in inseqs.fasta (-i) is chimeric by splitting into 
 4 (-n) roughly equal-sized fragments, and assigning taxonomy to each fragment
 using the blast taxon assigner (-t and -r, see qiime.assign_taxonomy.py -h).
 Chimeras are sequences where different fragments are assigned to different
 taxonomies. Seq ids for putative chimeras will be written to 
 inseqs_chimeric.txt (default, derived from -i) along with the taxonomies 
 assigned to each fragment.

   identify_chimeric_seqs.py -i inseqs.fasta -t id_to_taxonomy.txt -r refseqs.fasta -n 4
"""

required_options = [\
    make_option('-i', '--input_seqs_fp',
        help='Path to fasta file of sequences to be assigned')

]

chimera_detection_method_choices = ['blast_fragments']

optional_options = [\
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

def main():
    """Run chimera checker with given options>"""

    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

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
    input_seqs_fp = opts.input_seqs_fp
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
