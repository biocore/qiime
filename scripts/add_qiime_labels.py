#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.6.0"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Release"

from os.path import isdir 

from qiime.util import parse_command_line_parameters, get_options_lookup,\
 make_option, create_dir
from qiime.add_qiime_labels import add_qiime_labels


options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""Takes a directory and a mapping file of SampleIDs to fasta file names, combines all files that have valid fasta extensions into a single fasta file, with valid QIIME fasta labels."""
script_info['script_description']="""A tab separated text file with SampleIDs 
and fasta file names (just the file name itself, not the full or relative 
filepath) is used to generate a combined fasta file with valid
QIIME labels based upon the SampleIDs specified in the mapping file.

Example mapping file:
Sample.1	fasta_dir/seqs1.fna
Sample.2	fasta_dir/seqs2.fna

This script is to handle situations where fasta data comes already 
demultiplexed into a one fasta file per sample basis.  Apart from altering
the fasta label to add a QIIME compatible label at the beginning (example:
>FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_
could become 
>control.sample_1 FLP3FBN01ELBSX length=250 xy=1766_0111 region=1 run=R_2008_12_09_13_51_01_

Note that limited checking is done on the mapping file.  The only tests
are that every fasta file name is unique, and that SampleIDs are
MIMARKS compliant (alphanumeric and period characters only).  Duplicate 
SampleIDs are allowed, so care should be taken that there are no typos.

No changes are made to the sequences.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Specify fasta_dir as the input directory of fasta files, use the SampleID to fasta file mapping file example_mapping.txt, start enumerating with 1000000 following SampleIDs, and output the data to the directory combined_fasta""","""%prog -i fasta_dir -m example_mapping.txt -n 1000000 -o combined_fasta"""))
script_info['output_description']="""A combined_seqs.fasta file will be created in the output directory, with the sequences assigned to Sample.1 and Sample.2."""
script_info['required_options']= [\
    make_option('-m', '--mapping_fp',type='existing_filepath',
                help='SampleID to fasta file name mapping file filepath'),
    make_option('-i', '--fasta_dir',type='existing_dirpath',
                help='Directory of fasta files to combine and label.')
    
]
script_info['optional_options']= [\
    make_option('-o', '--output_dir',type='new_dirpath',
        help='Required output directory for log file and corrected mapping '+\
        'file, log file, and html file. [default: %default]', default="./"),
    make_option('-n', '--count_start',
        help='Specify the number to start enumerating sequence labels with. '+\
        '[default: %default]', default=0, type="int")
    ]
        
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
     parse_command_line_parameters(**script_info)
      
    mapping_fp = opts.mapping_fp
    fasta_dir = opts.fasta_dir
    output_dir = opts.output_dir
    count_start = int(opts.count_start)

    # Check input filepaths
    try:
        test_mapping_f = open(mapping_fp, "U")
    except IOError:
        raise IOError,("Cannot open mapping filepath "+\
         "%s, please check filepath and permissions." % mapping_fp)
         
    if not isdir(fasta_dir):
        raise IOError,("Specified fasta dir "+
         "%s, does not exist" % fasta_dir)
    
    # Create output directory, check path/access to mapping file
    create_dir(output_dir)
    
    add_qiime_labels(open(mapping_fp, "U"), fasta_dir, output_dir, count_start)


if __name__ == "__main__":
    main()
