#!/usr/bin/env python
# File created on 30 Aug 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from optparse import make_option
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,
                        merge_n_otu_tables)
from qiime.format import format_otu_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--input_fps',help='the otu tables'),\
 make_option('-o','--output_fp',help='the output otu table filepath'),\
]
script_info['optional_options'] = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
       
    input_fs = map(open,opts.input_fps.split(','))
    
    sample_names, otu_names, data, taxonomy = merge_n_otu_tables(input_fs)
    
    out_f = open(opts.output_fp,'w')
    out_f.write(format_otu_table(sample_names=sample_names,
                                 otu_names=otu_names,
                                 data=data,
                                 taxonomy=taxonomy))
    out_f.close()






if __name__ == "__main__":
    main()