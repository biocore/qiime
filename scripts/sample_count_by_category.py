#!/usr/bin/env python
# File created on 19 Dec 2012
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

from collections import defaultdict
from qiime.parse import parse_mapping_file
from qiime.util import parse_command_line_parameters, make_option
from qiime.sort import natsort

script_info = {}
script_info['brief_description'] = "Count the number of samples associated to a category value"
script_info['script_description'] = """Sum up the number of samples with each category value and print this information."""
script_info['script_usage'] = [("Example:","Count the number of samples associated with Treatment","""%prog -i $PWD/mapping.txt -c Treatment""")]
script_info['output_description']= """Two columns, the first being the category value and the second being the count. Output is to standard out. If there are unspecified values, the output category is identified as ***UNSPECIFIED***"""
script_info['required_options'] = [\
 make_option('-i','--input_fp',type="existing_filepath",help='the input filepath'),\
 make_option('-c','--category',type='string',help='the category to examine')
]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    map_data, header, comments = parse_mapping_file(opts.input_fp)
    
    if opts.category not in header:
        raise KeyError, "%s doesn't appear to exist in the mapping file!" % opts.category

    result = defaultdict(int)
    cat_idx = header.index(opts.category)
    for samp in map_data:
        result[samp[cat_idx]] += 1

    for cat_val in natsort(result):
        if not cat_val:
            print "***UNSPECIFIED***\t%d" % result[cat_val]
        else:
            print "%s\t%d" % (cat_val, result[cat_val])

if __name__ == "__main__":
    main()
