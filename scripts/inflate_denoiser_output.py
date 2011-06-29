#!/usr/bin/env python
# File created on 21 Apr 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from qiime.util import make_option
from itertools import chain
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import inflate_denoiser_output

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Inflate denoiser results so they can be passed directly to OTU pickers."
script_info['script_description'] = """Inflate denoiser results so they can be passed directly to pick_otus.py, parallel_pick_otus_uclust_ref.py, or pick_otus_through_otu_table.py. Note that the results of this script have not be abundance sorted, so they must be before being passed to the OTU picker. The uclust OTU pickers incorporate this abundance presorting by default.

The inflation process writes each centroid sequence n times, where n is the number of reads that cluster to that centroid, and writes each singleton once. Flowgram identifiers are mapped back to post-split_libraries identifiers in this process (i.e., identifiers in fasta fps).
"""
script_info['script_usage'] = [
 ("","Inflate the results of a single denoiser run.","inflate_denoiser_output.py -c centroids.fna -s singletons.fna -f seqs.fna -d denoiser_mapping.txt -o denoised_seqs.fna"),
 ("","Inflate the results of multiple denoiser runs to a single denoised_seqs.fna file.","inflate_denoiser_output.py -c centroids1.fna,centroids2.fna -s singletons1.fna,singletons2.fna -f seqs1.fna,seqs2.fna -d denoiser_mapping1.txt,denoiser_mapping2.txt -o denoised_seqs.fna")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-c','--centroid_fps',help='the centroid fasta filepaths'),
 make_option('-s','--singleton_fps',help='the singleton fasta filepaths'),
 make_option('-f','--fasta_fps',help='the input (to denoiser) fasta filepaths'),
 make_option('-d','--denoiser_map_fps',help='the denoiser map filepaths'),
 make_option('-o','--output_fasta_fp',help='the output fasta filepath'),
]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    centroid_seqs = \
     [MinimalFastaParser(open(e,'U')) for e in opts.centroid_fps.split(',')]
    singleton_seqs = \
     [MinimalFastaParser(open(e,'U')) for e in opts.singleton_fps.split(',')]
    fasta_seqs = \
     [MinimalFastaParser(open(e,'U')) for e in opts.fasta_fps.split(',')]
    denoiser_map_fs = \
     [open(e,'U') for e in opts.denoiser_map_fps.split(',')]
    output_fasta_fp = opts.output_fasta_fp
    
    output_f = open(opts.output_fasta_fp,'w')
    for s in inflate_denoiser_output(
                    chain.from_iterable(centroid_seqs),
                    chain.from_iterable(singleton_seqs),
                    chain.from_iterable(denoiser_map_fs),
                    chain.from_iterable(fasta_seqs)):
        output_f.write('>%s\n%s\n' % s)
    output_f.close()

if __name__ == "__main__":
    main()