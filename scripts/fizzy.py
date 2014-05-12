#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division
from qiime.util import parse_command_line_parameters, make_option
import qiime.fizzy as fizzy 

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison","Gail Rosen"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"
 

script_info = {}
script_info['brief_description'] = """Run feature selection on abundance data contained in a BIOM file."""
script_info['script_description'] ="""This script will run a feature selection algorithm on abundance data contained in a BIOM file given a mapping file. The current feature selection methods uses a forward search algorithm to select the features. The objective functions are based on information theory. At the moment, users are limited to the objective functions implemented in the PyFeast feature selection module. You can find a tutorial of fizzy at http://qiime.org/tutorials/feature_selection.html """
script_info['script_usage'] = [(
    """Run JMI feature selection on a BIOM file:""",
    """To perform feature selection the BIOM file, mapping file must be specified in advance, and the label column in the mapping file. Here we use JMI and select 15 features. """,
    """%prog -i data.biom -m map.txt -c Class -o result.txt -f JMI -k 15""")]
script_info['output_description']= """Text file containing the top features selected by the algorithm. """

script_info['required_options'] = [
    make_option('-c',
        '--column_label',
        type="string",
        help='column indicating the labels in the map file.'),
    make_option('-i', 
        '--input_fp',
        type="existing_filepath",
        help='input biom file'),
    make_option('-m',
        '--mapping_fp',
        type="existing_filepath",
        help='mapping file with labeling scheme'),
    make_option('-o',
        '--output_fp',
        type="new_filepath",
        help='the output file')
]
script_info['optional_options'] = [
    make_option('-k',
        '--n_select',
        type="int",
        help='number of features to select [default: %default]',
        default=15),
    make_option('-f',
        '--fs_method',
        type='choice', 
        help='feature selection method. valid options are ' 
        + ', '.join(fizzy.get_fs_methods()) + '. [default: %default]',
        choices=fizzy.get_fs_methods(),
        default='MIM')
]
script_info['version'] = __version__



def main():
    """
        main()
    """
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # run the fizzy feature selection routine
    selected_features = fizzy.run_feature_selection( 
        open(opts.input_fp,'U'), 
        open(opts.mapping_fp,'U'), 
        opts.column_label, 
        opts.fs_method, 
        opts.n_select)


    f = open(opts.output_fp, 'w')
    for sf in selected_features:
        f.write(sf + '\n')
    f.close()

if __name__ == "__main__":
    main()

