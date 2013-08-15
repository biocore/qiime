#!/usr/bin/env python
# File created on 15 Aug 2013
from __future__ import division

__author__ = "Will Van Treuren, Luke Ursell"
__copyright__ = "Copyright 2013s, The QIIME project"
__credits__ = ["Will Van Treuren, Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Will Van Treuren, Luke Ursell"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"


from qiime.util import parse_command_line_parameters, make_option
from qiime.ocs import sync_biom_and_mf


script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 #make_option('-i','--input_fp',type="existing_filepath",help='the input filepath'),\
]
script_info['optional_options'] = [\
 # Example optional option
 #make_option('-o','--output_dir',type="new_dirpath",help='the output directory [default: %default]'),\
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    # sync the mapping file and the biom file
    bt = parse_biom_table(open(opts.biom_fp))
    pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
    pmf, bt sync_biom_and_mf(pmf, bt)

    


if __name__ == "__main__":
    main()