#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena",
               "Catherine Lozupone", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

from os.path import exists, isdir, splitext, join
from os import makedirs, listdir

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option)
from qiime.principal_coordinates import pcoa

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Principal Coordinates Analysis (PCoA)"
script_info['script_description'] = ("Principal Coordinate Analysis (PCoA) is "
                                     "commonly used to compare groups of "
                                     "samples based on phylogenetic or "
                                     "count-based distance metrics (see "
                                     "section on beta_diversity.py).")
script_info['script_usage'] = [
    ("PCoA (Single File)",
     "For this script, the user supplies a distance matrix (i.e. resulting "
     "file from beta_diversity.py), along with the output filename (e.g.  "
     "beta_div_coords.txt), as follows:",
     "%prog -i beta_div.txt -o beta_div_coords.txt"),
    ("PCoA (Multiple Files):",
     "The script also functions in batch mode if a folder is supplied as input"
     " (e.g. from beta_diversity.py run in batch). No other files should be "
     "present in the input folder - only the distance matrix files to be "
     "analyzed. This script operates on every distance matrix file in the "
     "input directory and creates a corresponding principal coordinates "
     "results file in the output directory, e.g.:",
     "%prog -i beta_div_weighted_unifrac/ -o beta_div_weighted_pcoa_results/")
]
script_info['output_description'] = ("The resulting output file consists of "
                                     "the principal coordinate (PC) axes "
                                     "(columns) for each sample (rows). "
                                     "Pairs of PCs can then be graphed to view"
                                     " the relationships between samples. The "
                                     "bottom of the output file contains the "
                                     "eigenvalues and % variation explained "
                                     "for each PC. For more information of the"
                                     " file format, check the "
                                     "OrdinationResults class in the "
                                     "scikit-bio package "
                                     "(http://scikit-bio.org/)")
script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_path',
                help='path to the input distance matrix file(s) (i.e., the '
                     'output from beta_diversity.py). Is a directory for '
                     'batch processing and a filename for a single file '
                     'operation.'),
    make_option('-o', '--output_path', type='new_path',
                help='output path. directory for batch processing, filename '
                     'for single file operation'),
]

script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_path = opts.input_path
    output_path = opts.output_path

    if isdir(input_path):
        # Run PCoA on all distance matrices in the input dir
        # Create the output directory if it does not exists
        if not exists(output_path):
            makedirs(output_path)

        # Get all the filenames present in the input directory
        file_names = [fname for fname in listdir(input_path)
                      if not (fname.startswith('.') or isdir(fname))]

        # Loop through all the input files
        for fname in file_names:
            # Get the path to the input distance matrix
            infile = join(input_path, fname)

            # Run PCoA on the input distance matrix
            with open(infile, 'U') as lines:
                pcoa_scores = pcoa(lines)

            # Store the PCoA results on the output directory
            base_fname, ext = splitext(fname)
            out_file = join(output_path, 'pcoa_%s.txt' % base_fname)
            pcoa_scores.write(out_file)

    else:
        # Run PCoA on the input distance matrix
        with open(input_path, 'U') as f:
            pcoa_scores = pcoa(f)
        # Store the results in the output file
        pcoa_scores.write(output_path)


if __name__ == "__main__":
    main()
