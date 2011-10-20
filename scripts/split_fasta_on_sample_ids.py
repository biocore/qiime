#!/usr/bin/env python
# File created on 20 Oct 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from cogent.parse.fasta import MinimalFastaParser
from qiime.util import (parse_command_line_parameters, 
                        make_option, 
                        split_fasta_on_sample_ids_to_files)

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_fasta_fp',type="existing_filepath",help='the input fasta file to split'),
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: %default]'),\
]
script_info['optional_options'] = [\
 make_option('--buffer_size',type="int",default=500,
 help="the number of sequences to read into memory before writing to file (you usually won't need to change this) [default: %default]"),\
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    split_fasta_on_sample_ids_to_files(MinimalFastaParser(open(opts.input_fasta_fp,'U')),
                                       opts.output_dir,
                                       opts.buffer_size)


if __name__ == "__main__":
    main()