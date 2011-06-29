#!/usr/bin/env python
# File created on 14 Mar 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from qiime.util import make_option
from os.path import split
from qiime.util import parse_command_line_parameters, get_options_lookup, create_dir
from qiime.filter import split_otu_table_on_taxonomy

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Script to split a single OTU table into multiple tables based on the taxonomy at some user-specified depth."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Split seqs_otu_table.txt into taxon-specific OTU tables based on the third level in the taxonomy, and write the taxon-specific OTU tables to ./L3/","split_otu_table_by_taxonomy.py -i seqs_otu_table.txt -L 3 -o ./L3/")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--input_fp',help='the input otu table'),
 make_option('-L','--level',help='the level to split on',type="int"),
 make_option('-o','--output_dir',help='the output directory'),
]
script_info['optional_options'] = []
script_info['version'] = __version__

def split_otu_table_on_taxonomy_to_files(otu_table_fp,
                                         level,
                                         output_dir):
    """ Split OTU table by taxonomic level, writing otu tables to output dir
    """
    results = []
    create_dir(output_dir)
    filename_suffix = split(otu_table_fp)[1]
    for taxon_at_level, otu_table_str in\
     split_otu_table_on_taxonomy(open(otu_table_fp,'U'),level):
        taxon_fn = taxon_at_level.replace(';','-').replace(' ','_')
        output_fp = '%s/%s_%s' % (output_dir,taxon_fn,filename_suffix)
        output_f = open(output_fp,'w')
        output_f.write(otu_table_str)
        output_f.close()
        results.append((taxon_at_level,output_fp))
    return results

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    split_otu_table_on_taxonomy_to_files(opts.input_fp,
                                         opts.level,
                                         opts.output_dir)


if __name__ == "__main__":
    main()