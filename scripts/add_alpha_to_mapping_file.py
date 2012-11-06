#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Development"

from qiime.format import format_mapping_file
from qiime.parse import parse_matrix, parse_mapping_file
from qiime.util import parse_command_line_parameters, make_option
from qiime.add_alpha_to_mapping_file import add_alpha_diversity_values_to_mapping_file

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
make_option('-i','--alpha_fp',type="existing_filepath",\
    help='alpha_diversity.py output, with one or multiple metrics'),\
make_option('-m','--mapping_fp',type="existing_filepath",\
    help='input mapping file to modify'),\
]
script_info['optional_options'] = [\
make_option('-o','--output_mapping_fp',type="new_filepath",\
    help='output mapping file to modify [default: %default]',\
    default='mapping_file_with_alpha.txt'),\
make_option('-b','--number_of_bins',type="int",\
    help='Number of bins [default: %default]', default=4),\
make_option('-x','--missing_value_name',type="string",\
    help='Bin name for the sample identifiers that exist in the mapping file'+\
        ' (mapping_fp) but not in the alpha diversity file (alpha_fp)'+\
        '[default: %default]', default='N/A')\
]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    alpha_fp = opts.alpha_fp
    mapping_fp = opts.mapping_fp
    output_mapping_fp = opts.output_mapping_fp
    missing_value_name = opts.missing_value_name

    try:
        number_of_bins = int(opts.number_of_bins)
    except ValueError:
        raise ValueError, 'The number of bins must be an integer, not %s'\
            % opts.number_of_bins


    mapping_file_data, mapping_file_headers, comments = parse_mapping_file(\
        open(mapping_fp, 'U'))

    alpha_metrics, alpha_sample_ids, alpha_data = parse_matrix(\
        open(alpha_fp, 'U'))

    out_mapping_file_headers, out_mapping_file_data = \
        add_alpha_diversity_values_to_mapping_file(alpha_metrics,\
        alpha_sample_ids, alpha_data, mapping_file_headers,\
        mapping_file_data, number_of_bins, missing_value_name)

    lines = format_mapping_file(out_mapping_file_data, out_mapping_file_headers)

    fd_out = open(output_mapping_fp, 'w')
    fd_out.writelines(lines)
    fd_out.close()

if __name__ == "__main__":
    main()