#!/usr/bin/env python
# File created on 09 Feb 2010
#file filter_otus_by_sample.py

from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"
 
from optparse import make_option
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.filter_otus_by_sample import filter_samples,process_extract_samples
from qiime.make_3d_plots import create_dir
from cogent import LoadSeqs
from qiime.parse import fields_to_dict

options_lookup = get_options_lookup()

#filter_otus_by_sample.py
script_info={}
script_info['brief_description']="""Filter OTU mapping file and sequences by SampleIDs"""
script_info['script_description']="""This filter allows for the removal of sequences and OTUs containing user-specified Sample IDs, for instance, the removal of negative control samples. This script identifies OTUs containing the specified Sample IDs and removes its corresponding sequence from the sequence collection."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""The following command can be used, where all options are passed (using the resulting OTU file from pick_otus.py, FASTA file from split_libraries.py and removal of Sample1) with the resulting data being written to the output directory "filtered_otus/":""","""filter_otus_by_sample.py -i seqs_otus.txt -f seqs.fna -s Sample1 -o filtered_otus/"""))
script_info['output_description']="""As a result a new OTU and sequence file is generated and written to a randomly generated folder where the name of the folder starts with "filter_by_otus" Also included in the folder, is another FASTA file containing the removed sequences, leaving the user with 3 files."""

script_info['required_options']=[\
 options_lookup['input_fasta'],
 make_option('-i', '--input_otu_path', help='Path to OTU mapping file \
 containing sequence ids assigned to each OTU (i.e., resulting OTU file from \ pick_otus.py)'),
 make_option('-s', '--samples_to_extract', help='This is a list of sample \
ids, which should be removed from the OTU file')]

script_info['optional_options']=[\
options_lookup['output_dir']
]

script_info['version'] = __version__

def main():
    """opens files as necessary based on prefs"""
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    data = {}

    fasta_file = opts.input_fasta_fp

    # load the input alignment
    data['aln'] = LoadSeqs(fasta_file,aligned=False)

    #Load the otu file
    otu_path=opts.input_otu_path
    otu_f = open(otu_path, 'U')
    otus = fields_to_dict(otu_f)
    otu_f.close()

    data['otus']=otus
    #Determine which which samples to extract from representative seqs
    #and from otus file
    if opts.samples_to_extract:
      prefs=process_extract_samples(opts.samples_to_extract)

    filepath=opts.input_fasta_fp
    filename=filepath.strip().split('/')[-1]
    filename=filename.split('.')[0]

    dir_path = create_dir(opts.output_dir,'filtered_by_otus')

    try:
      action = filter_samples
    except NameError:
      action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
      action(prefs, data, dir_path,filename)

if __name__ == "__main__":
    main()