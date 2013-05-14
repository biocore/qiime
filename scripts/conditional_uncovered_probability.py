#!/usr/bin/env python
# File created on 1 April 2012
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Jose Antonio Navas Molina"]
#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Release" 

from qiime.util import make_option, parse_command_line_parameters
from qiime.alpha_diversity import (single_file_cup, 
    list_known_cup_metrics)
import os

#conditional_uncovered_probability.py
script_info={}

script_info['version'] = __version__
script_info['script_description'] = "Calculate the conditional uncovered\
 probability."

script_info['brief_description']="""Calculate the conditional uncovered\
 probability on each sample in an otu table."""
script_info['script_description']="""This script calculates the conditional\
 uncovered probability for each sample in an OTU table. It uses the methods\
 introduced in Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via\
 Poissonization: Accurate Measurements of the Microbial Unknown" PLoS 2011. 

Specifically, it computes a point estimate and a confidence interval using two \
different methods. Thus it can happen that the PE is actually outside of the CI. 

The CI method requires precomputed constants that depend on the lookahead, the\
 upper-to-lower bound ratio and the desired confidence.
We only provide these constants for some frequently used combinations. These\
 are (alpha:0.95, r=1..25)) for the the L and U interval types, and (alpha:0.9,\
 0.95, 0.99; f=10;  r=3..25,30,40,50). Also, there are a few hand picked\
 special cases:

f=2 and r=50 and alpha=0.95
f=2 and r=33 and alpha=0.95
f=1.5 and r=100 and alpha=0.95
f=1.5 and r=94 and alpha=0.95
f=2.5 and r=19 and alpha=0.95

"""

script_info['script_usage']=[]

script_info['script_usage'].append(("""Default case:""",
"""To calculate the cond. uncovered probability with the default values,\
 you can use the following command: """,
"""%prog -i otu_table.biom -o cup.txt"""))

script_info['script_usage'].append(("""Change lookahead:""",
"""To change the accuracy of the prediction change the lookahead value. Larger\
 values of r lead to more precise predictions, but might be unfeasable for\
 small samples. For deeply sequenced samples, try increasing r to 50: """,
"""%prog -i otu_table.biom -o cup_r50.txt -r 50"""))

script_info['script_usage'].append(("""Change the interval type:""",
"""To change the confidence interval type to a lower bound prediction, while\
 the upper bound is set to 1 use: """,
"""%prog -i otu_table.biom -o cup_lower_bound.txt -c L"""))

script_info['output_description']="""The resulting file(s) is a tab-delimited\
 text file, where the columns correspond to estimates of the cond. uncovered\
 probability and the rows correspond to samples. The output file is compatible\
 with the alpha_diversity output files and thus could be tied into the\
s rarefaction workflow.

Example Output:

====== ======= ============= ================
\      PE      Lower Bound   Upper Bound
====== ======= ============= ================
PC.354 0.111   0.0245        0.245
PC.124 0.001   0.000564      0.00564
====== ======= ============= ================

"""

script_info['required_options']=[]

CI_TYPES=['ULCL', 'ULCU', 'U', 'L']
script_info['optional_options']=[\
    make_option('-i', '--input_path',
                help='Input OTU table filepath. [default: %default]',
                type='existing_path'),
    make_option('-o', '--output_path',
             help='Output filepath to store the predictions. [default: %default]',
             type='new_path'),
    make_option('-r', '--look_ahead',
                help='Number of unobserved, new colors necessary for prediction.'
                ' [default: %default]', default=25,
                type='int'),
    make_option('-c', '--ci_type',
                help='Type of confidence interval.  Choice of '+\
                    ", ".join(CI_TYPES) +' [default: %default]', default='ULCL',
                type='choice', choices=CI_TYPES),
    make_option('-a', '--alpha',
                help='Desired confidence level for CI prediction.'+
                ' [default: %default]', default=0.95,
                type='float'),
    make_option('-f', '--f_ratio',
                help='Upper to lower bound ratio for CI prediction.' +
                ' [default: %default]', default=10.0,
                type='float'),
    make_option('-m', '--metrics', default='lladser_pe,lladser_ci',
                type='multiple_choice', mchoices=list_known_cup_metrics(),
                help='CUP metric(s) to use. A comma-separated list should' +\
                   ' be provided when multiple metrics are specified. [default: %default]'), 
    make_option('-s', '--show_metrics', action='store_true', 
                dest="show_metrics",
                help='Show the available CUP metrics and exit.'),
]

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.show_metrics:
        print("Known metrics are: %s\n" \
              % (', '.join(list_known_cup_metrics()),))
        exit(0)

    almost_required_options = ['input_path', 'output_path']
    for option in almost_required_options:
        if getattr(opts,option) == None:
            option_parser.error('Required option --%s omitted.' % option)
    
    if os.path.isfile(opts.input_path):
      try:
          f = open(opts.output_path, 'w')
      except IOError:
          exit("ioerror, couldn't create output file") 
      f.close()
      
      single_file_cup(opts.input_path, opts.metrics, opts.output_path,
                      opts.look_ahead, opts.alpha, opts.f_ratio, opts.ci_type)

    else:
      exit("io error, input path not valid. does it exist?")

if __name__ == "__main__":
    main()
