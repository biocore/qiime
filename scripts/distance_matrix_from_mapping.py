#!/usr/bin/env python
# File created on 23 Jul 2011
# based on beta_diversity.py
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena", "Andrew J. King"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"


from qiime.util import parse_command_line_parameters, make_option
from qiime.format import format_distance_matrix
from qiime.distance_matrix_from_mapping import compute_distance_matrix_from_metadata, calculate_dist_vincenty
from qiime.parse import parse_mapping_file_to_dict
import os.path

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')


script_info = {}
script_info[
    'brief_description'] = """Calculate the pairwise dissimilarity on one column of a mappping file"""
script_info['script_description'] = """The input for this script is a mapping file and the name of a column, it has to be numeric, from which a distance matrix will be created. The output of this script is a distance matrix containing a dissimilarity value for each pairwise comparison.

As this is a univariate procedure only one metric is supported: d = c-b."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Pairwise dissimilarity:",
     "To calculate the distance matrix (using euclidean distance) on a column of the mapping file, where the results are output to DOB.txt, use the following command:",
     "%prog -i Fasting_Map.txt -c DOB"))
script_info['script_usage'].append(
    ("Pairwise dissimilarity using the Vincenty formula for distance between two Latitude/Longitude points:",
     "To calculate the distance matrix (using Vincenty formula) on a column of the mapping file, where the results are output to lat_long.txt, use the following command:",
     "%prog -i lat_long.txt -c Latitute,Longitude -o lat_long_dtx_matrix.txt"))
script_info[
    'output_description'] = """The output of distance_matrix_from_mapping.py is a file containing a distance matrix between rows corresponding to a pair of columns in a mapping file."""
script_info['required_options'] = [
    make_option('-i', '--input_path',
                help='Mapping filepath.', type='existing_path'),
    make_option('-c', '--column', type='string',
                help="string containing the name of the column in the mapping file, e.g. 'DOB'. If you pass two colums separated by a comma (e.g. 'Latitude,Longitud') the script will calculate the Vincenty formula (WGS-84) for distance between two Latitude/Longitude points."),
]
script_info['optional_options'] = [
    make_option('-o', '--output_fp',
                help="Output path to store the distance matrix. [default=%default]",
                type='new_filepath', default="map_distance_matrix.txt")
]
script_info['option_label'] = {'input_path': 'Mapping filepath',
                               'column': 'List of samples for compute',
                               'output_fp': 'Output filename'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    data, comments = parse_mapping_file_to_dict(opts.input_path)
    column_headers = []
    if ',' not in opts.column:
        column_data = []
        column_name = opts.column
        for i in data:
            if column_name not in data[i]:
                raise ValueError(
                    "No column: '%s' in the mapping file. Existing columns are: %s" %
                    (column_name, data[i].keys()))

            try:
                column_data.append(float(data[i][opts.column]))
            except ValueError:
                raise ValueError(
                    "All the values in the column '%s' must be numeric but '%s' has '%s'" %
                    (column_name, i, data[i][column_name]))

            column_headers.append(i)
        dtx_mtx = compute_distance_matrix_from_metadata(column_data)
    else:
        latitudes = []
        longitudes = []
        try:
            latitude, longitude = opts.column.split(',')
        except ValueError:
            raise ValueError(
                "This script accepts a maximum of 2 colums separated by comma and you passed: %s" %
                (opts.column))

        for i in data:
            if latitude not in data[i] or longitude not in data[i]:
                raise ValueError(
                    "One of these columns or both do not exist: '%s' or '%s' in the mapping file. Existing columns are: %s" %
                    (latitude, longitude, data[i].keys()))

            try:
                latitudes.append(float(data[i][latitude]))
                longitudes.append(float(data[i][longitude]))
            except ValueError:
                raise ValueError(
                    "All the values in the columnd '%s' & '%s' must be numeric but '%s' has '%s'" %
                    (latitude, longitude, i, data[i][column_name]))

            column_headers.append(i)

        dtx_mtx = calculate_dist_vincenty(latitudes, longitudes)

    dtx_txt = format_distance_matrix(column_headers, dtx_mtx)

    outfilepath = os.path.join(opts.output_fp)
    f = open(outfilepath, 'w')
    f.write(dtx_txt)
    f.close()


if __name__ == "__main__":
    main()
