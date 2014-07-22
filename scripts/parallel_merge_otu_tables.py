#!/usr/bin/env python

from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"


from qiime.util import (parse_command_line_parameters, load_qiime_config,
                        get_options_lookup, make_option)
from qiime.parallel.merge_otus import ParallelMergeOtus

qiime_config = load_qiime_config()
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Parallel merge BIOM tables"
script_info['script_description'] = (
    "This script works like the merge_otu_tables.py script, but is intended "
    "to make use of multicore/multiprocessor environments to perform analyses"
    " in parallel.")
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Example",
     "Merge the OTU tables $PWD/t1.biom,$PWD/t2.biom,$PWD/t3.biom,$PWD/t4.biom"
     " and write the resulting output table to the $PWD/merged/ directory.",
     "%prog -i $PWD/t1.biom,$PWD/t2.biom,$PWD/t3.biom,$PWD/t4.biom -o "
     "$PWD/merged/"))
script_info['output_description'] = (
    "The output consists of many files (i.e. merged_table.biom, "
    "merged_table.log and all intermediate merge tables). The .biom file "
    "contains the result of merging the individual BIOM tables. The resulting "
    ".log file contains a list of parameters passed to this script along with "
    "the output location of the resulting .txt file, the dependency hierarchy "
    "and runtime information for each individual merge.")

script_info['required_options'] = [
    make_option('-i', '--input_fps', type='existing_filepaths',
                help='the otu tables in biom format (comma-separated)'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output otu table directory path')]

script_info['optional_options'] = [
    options_lookup['retain_temp_files'],
    options_lookup['suppress_submit_jobs'],
    options_lookup['suppress_blocking']]

script_info['version'] = __version__


if __name__ == '__main__':
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    parallel_runner = ParallelMergeOtus(
        retain_temp_files=opts.retain_temp_files,
        block=not opts.suppress_blocking)
    parallel_runner(opts.input_fps, opts.output_dir)
