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
from os import makedirs
from os.path import exists, splitext, split, isdir

from qiime.align_seqs import alignment_module_names,alignment_method_constructors,\
    pairwise_alignment_methods, CogentAligner
from qiime.util import load_qiime_config

script_description = """Align sequences using a variety of alignment methods. """

script_usage = """
Align 10_seq.fasta (-i) using muscle (-m). Output files will be stored in 
./muscle_aligned/ 
  align_seqs.py -i 10_seq.fasta -m muscle
 
Align 10_seq.fasta (-i) with pynast (default) against template_aln.fasta (-t).
Output files will be stored in ./pynast_aligned/
  align_seqs.py -i 10_seq.fasta -t template_aln.fasta

Align 10_seq.fasta (-i) with infernal against template_aln.sto (-t).
Output files will be stored in ./infernal_aligned/
  align_seqs.py -i 10_seq.fasta -t template_aln.sto -m infernal
"""

qiime_config = load_qiime_config()

required_options = [\
   make_option('-i','--input_fasta_fp',\
        help='Input sequences to align')
]


if qiime_config['pynast_template_alignment_fp']:
    template_fp_default_help = '[default: %default]'
else:
    template_fp_default_help = '[REQUIRED if -m pynast or -m infernal]'

blast_db_default_help =\
     qiime_config['pynast_template_alignment_blastdb'] or\
      'created on-the-fly from template_alignment'

alignment_method_choices = \
    alignment_method_constructors.keys() + alignment_module_names.keys()
pairwise_alignment_method_choices = pairwise_alignment_methods.keys()

optional_options = [\
    make_option('-t','--template_fp',\
          type='string',dest='template_fp',help='Filepath for '+\
          'template against %s' % template_fp_default_help,\
          default=qiime_config['pynast_template_alignment_fp']),

    make_option('-m','--alignment_method',\
          type='choice',help='Method for aligning'+\
          ' sequences. Valid choices are: ' +\
          ', '.join(alignment_method_choices) + ' [default: %default]',
          choices=alignment_method_choices,\
          default='pynast'),
          
    make_option('-a','--pairwise_alignment_method',\
          type='choice',help='method for performing pairwise ' +\
          'alignment in PyNAST. Valid choices are '+\
          ', '.join(pairwise_alignment_method_choices) +\
          ' [default: %default]',\
          choices=pairwise_alignment_method_choices,\
          default='blast'),

    make_option('-d','--blast_db',\
          dest='blast_db',help='Database to blast against when -m pynast '+\
          '[default: %s]' % blast_db_default_help,\
          default=qiime_config['pynast_template_alignment_blastdb']),
          
    make_option('-o','--output_dir',\
          help='Path to store '+\
          'result file [default: <ALIGNMENT_METHOD>_aligned]'),
          
    make_option('-e','--min_length',\
          type='int',help='Minimum sequence '+\
          'length to include in alignment [default: %default]',\
           default=150),
          
    make_option('-p','--min_percent_id',\
          type='float',help='Minimum percent '+\
          'sequence identity to closest blast hit to include sequence in'+\
          ' alignment [default: %default]', default=0.75)
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    
    if opts.min_percent_id <= 1.0:
        opts.min_percent_id *= 100
        
    if not (1.0 <= opts.min_percent_id <= 100.0):
        option_parser.error('Minimum percent sequence identity must be' +\
        ' between 1.0 and 100.0: %2.2f' % opts.min_percent_id)
        
    if not opts.template_fp and opts.alignment_method == 'pynast':
        option_parser.error('PyNAST requires a template alignment to be passed via -t')

    input_seqs_filepath = opts.input_fasta_fp
    alignment_method = opts.alignment_method
    output_dir = opts.output_dir or alignment_method + '_aligned'
    #TODO: should check if output dir is writable before doing actual computation

    if (exists(output_dir)):
        if isdir(output_dir):
            #Looks good, dir is there
            # maybe issue a warning to log file?
            # should check if dir is writable
            pass
        else:
            #File with same name
            raise OSError,"File with same name as outdir exists: %s" % output_dir
    else:
        #no dir there, make it
        try:
            makedirs(output_dir)
        except OSError:
            #re-raise error, but slightly more informative 
            raise OSError,"Could not create output directory %s" % output_dir
    
    fpath, ext = splitext(input_seqs_filepath)
    input_dir, fname = split(fpath)
    
    result_path = output_dir + '/' + fname + "_aligned" + ext
    log_path = output_dir + '/' + fname + "_log.txt"
    failure_path = output_dir + '/' + fname + "_failures.fasta"
 
    if alignment_method in alignment_method_constructors:
        # try/except was causing problems here, so replacing with
        # an explicit check
        # define the aligner params
        aligner_constructor =\
         alignment_method_constructors[alignment_method]
        aligner_type = alignment_method
        params = {'min_len':opts.min_length,\
                  'min_pct':opts.min_percent_id,\
                  'template_filepath':opts.template_fp,\
                  'blast_db':opts.blast_db,
                  'pairwise_alignment_method': opts.pairwise_alignment_method}
        # build the aligner object
        aligner = aligner_constructor(params)
        # apply the aligner
        aligner(input_seqs_filepath,result_path=result_path,\
         log_path=log_path,failure_path=failure_path)
    else:
        # define the aligner params
        aligner = CogentAligner({
        'Module':alignment_module_names.get(alignment_method, 'Unknown'),
        'Method':alignment_method
        })
        # build the aligner
        aligner_type = 'Cogent'
        # apply the aligner
        aligner(result_path, seq_path=input_seqs_filepath,
            log_path=log_path)

if __name__ == "__main__":
    main()
