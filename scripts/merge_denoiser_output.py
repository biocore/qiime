#!/usr/bin/env python
# File created on 01 Apr 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Jens Reeder"
__email__ = "jreeder@colorado.edu"
__status__ = "Release"
 
from qiime.util import make_option
from cogent.app.util import ApplicationNotFoundError
from qiime.util import parse_command_line_parameters, get_options_lookup, create_dir

# Adapted from align_seqs.py
# Load Denoiser if it's available. If it's not, skip it if not but set up
# to raise errors if the user tries to use it.
try:
    from Denoiser.denoise_postprocess import post_process

except ImportError:
    def raise_denoiser_not_found_error(*args, **kwargs):
        raise ApplicationNotFoundError,\
         "Denoiser cannot be found.\nIs it installed? Is it in your $PYTHONPATH?"+\
         "\nYou can obtain the Denoiser from http://www.microbio.me/denoiser .\n"
    # set functions which cannot be imported to raise_denoiser_not_found_error
    post_process= raise_denoiser_not_found_error

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Merge the output of denoising step back into QIIME"
script_info['script_description'] = """
This script combines the output of the denoising step with the OTU picker results.
This script has to be run after denoising and OTU picking to combine the denoiser clusters and OTU clusters.
The input to the script is:

   DENOISER_MAP_FILE: the cluster mapping from the denoiser, usually denoiser_mapping.txt in the denoiser output directory
   
   DENOISED_FASTA_FP: path to the output FASTA file of the denoiser (centroids.fasta and singletons.fasta in the denoiser output dir. Concatenate if you want to include both)

   OTU_PICKER_MAP_FILE: path to OTU mapping file from OTU picking on the

   FASTA_FP:  path to FASTA input file, thus has to be the output of split_libraries.py
             

See the tutorial on 454 Denoising in the QIIME tutorials on how to use this script properly.
"""
script_info['script_usage'] = [("Example usage",
                                "Merge the output of denoising (denoised_seqs.fasta and denoiser_mapping.txt) the OTU picker results on denoised_seqs.fasta (uclust_picked_otus/denoised_seqs_otus.txt) and replace the read IDs with the sampleIDs from the output of split_libraries.py (seqs.fna)\n\n",
                                "merge_denoiser_output.py -f seqs.fna -d denoised_seqs.fasta -m denoiser_mapping.txt -p uclust_picked_otus/denoised_seqs_otus.txt")]
script_info['output_description']= """
Two files are generated:

1. A new otu mapping file with the combined clusters from denoising and OTU picking,

2. a FASTA file with one sequence per OTU.

These two files need to be provided to pick_rep_set using -m first"""

script_info['required_options'] = [\
    
    make_option('-m','--map_file',action='store',\
                    type='string',dest='denoiser_map_file',\
                    help='path to denoiser cluster mapping file '+\
                    '[REQUIRED]'),

    make_option('-p','--otu_picker_map_file', action='store',\
                    type='string',dest='otu_picker_map_file',\
                    help='path to OTU mapping file from OTU picker '+\
                    '[REQUIRED]'),
    
    make_option('-f','--fasta_fp',action='store',\
                    type='string',dest='fasta_fp',help='path to FASTA file, '+\
                    'output of split_libraries.py'+\
                    ' [REQUIRED]'),

    make_option('-d','--denoised_fasta_fp',action='store',\
                    type='string',dest='denoised_fasta_fp',\
                    help='path to denoised fasta file '+\
                    '[REQUIRED]')
]

script_info['optional_options'] = [\

    make_option('-o','--output_dir',action='store',\
                    default="Denoiser_out_otu_picked/",\
                    type='string',dest='output_dir',help='path to output'+\
                    ' directory [default: %default]')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    create_dir(opts.output_dir, fail_on_exist=False)
           
    post_process(opts.fasta_fp, opts.denoiser_map_file, opts.denoised_fasta_fp,
                 opts.otu_picker_map_file, opts.output_dir)
    
if __name__ == "__main__":
    main()
