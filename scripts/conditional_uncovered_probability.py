#!/usr/bin/env python
# File created on 1 April 2012
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Jose Antonio Navas Molina", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"

from qiime.util import make_option, parse_command_line_parameters
from qiime.alpha_diversity import (single_file_cup,
                                   list_known_cup_metrics)
import os

# conditional_uncovered_probability.py
script_info = {}

script_info['version'] = __version__
script_info['script_description'] = "Calculate the conditional uncovered\
 probability."

script_info['brief_description'] = """Calculate the conditional uncovered\
 probability on each sample in an otu table."""
script_info['script_description'] = \
"""This script calculates the conditional uncovered probability for each sample
in an OTU table. It uses the methods introduced in Lladser, Gouet, and Reeder,
"Extrapolation of Urn Models via Poissonization: Accurate Measurements of the
Microbial Unknown" PLoS 2011.

Specifically, it computes a point estimate and a confidence interval using two
different methods. Thus it can happen that the PE is actually outside of the
CI.

We only provide the ability to generate 95% (alpha=0.95) CIs. The CIs are ULCL
CIs; they provide an upper and lower bound, where the lower bound is
conservative. The CIs are constructed using an upper-to-lower bound ratio of
10.

The CI method requires precomputed constants that depend on the lookahead. We
only provide constants for r=3..25,30,40,50.

"""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Default case:",
     "To calculate the cond. uncovered probability with the default values, "
     "you can use the following command:",
     "%prog -i otu_table.biom -o cup.txt"))

script_info['script_usage'].append(
    ("Change lookahead:",
     "To change the accuracy of the prediction change the lookahead value. "
     "Larger values of r lead to more precise predictions, but might be "
     "unfeasable for small samples. For deeply sequenced samples, try "
     "increasing r to 50:",
     "%prog -i otu_table.biom -o cup_r50.txt -r 50"))

script_info['output_description'] = \
"""The resulting file(s) is a tab-delimited text file, where the columns
correspond to estimates of the cond. uncovered probability and the rows
correspond to samples. The output file is compatible with the alpha_diversity
output files and thus could be tied into the rarefaction workflow.

Example Output:

====== ======= ============= ================
\      PE      Lower Bound   Upper Bound
====== ======= ============= ================
PC.354 0.111   0.0245        0.245
PC.124 0.001   0.000564      0.00564
====== ======= ============= ================

"""

script_info['required_options'] = []

script_info['optional_options'] = [
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
    make_option('-m', '--metrics', default='lladser_pe,lladser_ci',
                type='multiple_choice', mchoices=list_known_cup_metrics(),
                help='CUP metric(s) to use. A comma-separated list should' +
                ' be provided when multiple metrics are specified. [default: %default]'),
    make_option('-s', '--show_metrics', action='store_true',
                dest="show_metrics",
                help='Show the available CUP metrics and exit.')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.show_metrics:
        print("Known metrics are: %s\n"
              % (', '.join(list_known_cup_metrics()),))
        exit(0)

    almost_required_options = ['input_path', 'output_path']
    for option in almost_required_options:
        if getattr(opts, option) is None:
            option_parser.error('Required option --%s omitted.' % option)

    if os.path.isfile(opts.input_path):
        try:
            f = open(opts.output_path, 'w')
        except IOError:
            exit("ioerror, couldn't create output file")
        f.close()

        single_file_cup(opts.input_path, opts.metrics, opts.output_path,
                        opts.look_ahead)

    else:
        exit("io error, input path not valid. does it exist?")

if __name__ == "__main__":
    main()
