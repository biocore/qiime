 #!/usr/bin/env python 

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Jens Reeder", "Rob Knight"]#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Release"

from os.path import exists

from cogent.util.misc import create_dir

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.denoiser.denoise_postprocess import combine_mappings

options_lookup = get_options_lookup()

#denoiser_postprocess.py
script_info={}
script_info['brief_description']="""Merge denoiser output with OTU picker output"""
script_info['script_description']="""Combine mapping files from denoiser and OTU picker, such that we have a combined mapping file that can be used for the subsequent steps of Qiime. Also replace flowgram identifiers with IDs assigned by split_libraries.py."""

script_info['script_usage'] = [ \
    ( "",
      """Combine denoiser output with output of QIIME OTU picker, put results into Outdir:""",
      """%prog -f seqs.fna -d denoised.fasta -m denoiser_mapping.txt -p cdhit_picked_otus/denoised_otus.txt -v -o Outdir""" )
    ]

script_info['output_description']=""" The output of denoiser_postprocess.py are two files:

denoised_otu_map.txt: In this mapping the read/flowgram IDs are
                replaced by sample_id from the split_libraries.py
                fasta file. Also, the lists for each OTU are sorted
                such that the largest cluster from denosing appears
                first. This will be important for the next step,
                picking representative sequences.

denoised_all.fasta: A fasta sequence where the header lines are
                updated with the sample_ids as assigned by split_libraries.py.
"""

script_info['required_options']=[\

    make_option('-p','--otu_picker_map_file',action='store',\
                    type='string',dest='otu_picker_map_file',\
                    help='path to OTU picker mapping file '+\
                    '[REQUIRED]', default= None),

    make_option('-f','--fasta_fp',action='store',\
                    type='string',dest='fasta_fp',\
                    help='path to fasta input file, '+\
                    'output of split_libraries.py'+\
                    ' [REQUIRED]', default=None),

    make_option('-d','--denoised_fasta_fp',action='store',\
                    type='string',dest='denoised_fasta_fp',\
                    help='path to denoised fasta file '+\
                    '[REQUIRED]', default=None)
    ]

script_info['optional_options']=[ \
    
    make_option('-o','--output_dir',action='store',\
                    type='string',dest='output_dir',\
                    help='path to output'+\
                    ' directory [default: %default]',\
                    default= "Denoiser_out_otu_picked/"),

    make_option('-m','--map_file',action='store',\
                    type='string',dest='denoiser_map_file',\
                    help='path to denoiser mapping file '+\
                    '[default: %default]', default="denoiser_mapping.txt")
    ]

script_info['version'] = __version__

def main(commandline_args=None):
    option_parser, opts, args = parse_command_line_parameters(**script_info)
 
    #check for missing files
    required_files = [opts.denoiser_map_file, opts.otu_picker_map_file,
                      opts.fasta_fp, opts.denoised_fasta_fp]
    if (not all(required_files) or not all(map(exists, required_files))):
        option_parser.error('Missing input files.')

    create_dir(opts.output_dir, fail_on_exist=False)
           
    combine_mappings(open(opts.fasta_fp), open(opts.denoiser_map_file),
                     open(opts.denoised_fasta_fp),
                     open(opts.otu_picker_map_file), opts.output_dir)

if __name__ == "__main__":
    main()
