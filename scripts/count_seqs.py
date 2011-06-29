#!/usr/bin/env python
# File created on 29 May 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from glob import glob
from qiime.util import (parse_command_line_parameters, 
                        make_option, 
                        count_seqs_in_filepaths)

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [\
 ("",
  "Count the sequences in a fasta file and write results to stdout.",
  "%prog -i in.fasta"),
 ("",
  "Count the sequences in two fasta files and write results to file.",
  "%prog -i in1.fasta,in2.fasta -o seq_counts.txt"),
  ("",
   "Count the sequences all .fasta files in current directory and write results to stdout. Note that -i option must be quoted.",
   "%prog -i \"*.fasta\"")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--input_fps',
        help='the input filepaths (comma-separated)'),
]
script_info['optional_options'] = [
 # Example optional option
 make_option('-o','--output_fp',type="new_filepath",
  help='the output filepath [default: write to stdout]'),\
 make_option('--suppress_errors',action='store_true',\
        help='Suppress warnings about missing files [default: %default]',
        default=False)
]
script_info['version'] = __version__

def format_output(count_data, total, inaccessible_filepaths, suppress_errors=False):
    """ Output formatter """
    lines = ['']
    for c in count_data:
        lines.append('%d  : %s (Sequence lengths (mean +/- std): %1.4f +/- %1.4f)' % 
                     (c[0][0],c[1],c[0][1],c[0][2]))
    lines.append('%d  : Total' % total)
    
    if inaccessible_filepaths and not suppress_errors:
        lines.append('')
        lines.append(\
         'Some files were not accessible. Do they exist? Do you have read permission?')
        for inaccessible_filepath in inaccessible_filepaths:
            lines.append('  %s' % inaccessible_filepath)
        lines.append('')
    return '\n'.join(lines)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    suppress_errors = opts.suppress_errors
    input_fps = []
    for e in opts.input_fps.split(','):
        input_fps.extend(glob(e))
    input_fps = set(input_fps)
    output_fp = opts.output_fp

    count_data, total, inaccessible_filepaths = count_seqs_in_filepaths(input_fps)
    r = format_output(count_data, total, inaccessible_filepaths, suppress_errors)
    
    if opts.output_fp:
        f = open(output_fp,'w')
        f.write(r)
        f.close()
    else:
        print r


if __name__ == "__main__":
    main()