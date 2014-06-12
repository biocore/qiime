#!/usr/bin/env python
# File created on 09 Feb 2010
# file make_otu_heatmap_html.py

from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jesse Stombaugh", "Jose Carlos Clemente Litran",
               "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

import os
import shutil
from sys import exit

from numpy import maximum

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.make_otu_heatmap_html import (generate_heatmap_plots,
                                         get_otu_counts)
from qiime.make_otu_heatmap import get_log_transform
from qiime.util import get_qiime_project_dir, create_dir
from qiime.parse import parse_mapping_file
from qiime.parse import parse_newick, PhyloNode

options_lookup = get_options_lookup()

# make_otu_heatmap_html.py
script_info = {}
script_info['brief_description'] = "Make heatmap of OTU table"
script_info['script_description'] = (
    "Create an interactive OTU heatmap from an OTU table. This script parses "
    "the OTU count table and filters the table by counts per otu "
    "(user-specified), then converts the table into a javascript array, which "
    "can be loaded into a web application. The OTU heatmap displays raw OTU "
    "counts per sample, where the counts are colored based on the "
    "contribution of each OTU to the total OTU count present in that sample "
    "(blue: contributes low percentage of OTUs to sample; red: contributes "
    "high percentage of OTUs). This web application allows the user to filter "
    "the otu table by number of counts per otu. The user also has the ability "
    "to view the table based on taxonomy assignment. Additional features "
    "include: the ability to drag rows (up and down) by clicking and dragging "
    "on the row headers; and the ability to zoom in on parts of the heatmap "
    "by clicking on the counts within the heatmap.")
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Generate an OTU heatmap",
     "By using the default values (-n 5), you can then use the code as "
     "follows:",
     "%prog -i otu_table.biom -o heatmap/"))

script_info['script_usage'].append(
    ("Generate a filtered OTU heatmap",
     "If you would like to filter the OTU table by a different number of "
     "counts per OTU (i.e., 10):",
     "%prog -i otu_table.biom -n 10 -o heatmap_mc10/"))

script_info['script_usage'].append(
    ("Generate a sample-sorted OTU heatmap",
     "If you would like to sort the heatmap by Sample IDs then you should "
     "supply the mapping file, as follows:",
     "%prog -i otu_table.biom -o heatmap_sample_sorted -m Fasting_Map.txt"))

script_info['script_usage'].append(
    ("Generate a sample and OTU-sorted OTU heatmap",
     "If you would like to sort the heatmap by Sample IDs and the tips in the "
     "tree, you can supply a tree as follows:",
     "%prog -i otu_table.biom -o heatmap_sample_otu_sorted -m Fasting_Map.txt "
     "-t rep_set.tre"))

script_info['output_description'] = (
    "The interactive heatmap is located in OUTPUT_DIR/otu_table.html where "
    "OUTPUT_DIR is specified as -o. Safari is recommended for viewing the "
    "OTU Heatmap, since the HTML table generation is much faster than "
    "Firefox (as of this writing).")

script_info['required_options'] = [
    options_lookup['otu_table_as_primary_input'],
    make_option('-o', '--output_dir', type="new_dirpath",
                help='path to the output directory'),
]
script_info['optional_options'] = [
    make_option('-n', '--num_otu_hits',
                help='Only include OTUs with at least this many sequences. '
                '[default: %default]', default=5, type='int'),
    make_option('-t', '--tree', help='Path to newick tree where OTUs are '
                'tips, used for sorting OTUs in the heatmap', default=None,
                type='existing_filepath'),
    make_option('-m', '--map_fname',
                help='Input metadata mapping filepath, used for sorting '
                'samples in the heatmap', default=None,
                type='existing_filepath'),
    make_option('--sample_tree', help='Path to newick tree where samples are '
                'tips (e.g, output from upgma_cluster.py) used for sorting '
                'samples in the heatmap. If both this and the metadata '
                'mapping file are provided, the mapping file will be ignored.',
                default=None, type='existing_filepath'),
    make_option('--log_transform', action="store_true", help='Log-transform '
                'the data. All zeros will be set to a small value (default is '
                '1/2 of the smallest non-zero entry). Data will be translated '
                'to be non-negative after log transform and the num_otu_hits '
                'will be set to 0.', default=False)
]
script_info['option_label'] = {
    'otu_table_as_primary_input': 'OTU table filepath',
    'map_fname': 'QIIME-formatted mapping filepath',
    'num_otu_hits': '# of sequences',
    'output_dir': 'Output directory',
    'tree': 'Newick tree filepath',
    'sample_tree': 'Newick tree filepath (containing samples on tips)',
    'log_transform': 'Perform log transformation'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    data = {}

    # Open and get coord data
    otu_table = get_otu_counts(opts.otu_table_fp)
    # determine whether fractional values are present in OTU table
    num_otu_hits = opts.num_otu_hits
    if opts.log_transform:
        otu_table = get_log_transform(otu_table)
        num_otu_hits = 0

    fractional_values = False
    max_val = -1
    for val in otu_table.iter_data(axis='observation'):
        max_val = maximum(max_val, val.max())

    # the data cannot be of mixed types: if one is float, all are float
    fractional_values = (
        max_val.dtype.name == 'float32' or max_val.dtype.name == 'float64')

    if fractional_values and max_val <= 1:
        if num_otu_hits > 0:
            print ("Warning: OTU table appears to be using relative "
                   "abundances and num_otu_hits was set to %d. Setting "
                   "num_otu_hits to 0." % num_otu_hits)
            num_otu_hits = 0

    filepath = opts.otu_table_fp
    filename = filepath.strip().split('/')[-1].split('.')[0]

    dir_path = opts.output_dir
    create_dir(dir_path)

    js_dir_path = os.path.join(dir_path, 'js')
    create_dir(js_dir_path)

    qiime_dir = get_qiime_project_dir()

    js_path = os.path.join(qiime_dir, 'qiime/support_files/js')

    shutil.copyfile(os.path.join(js_path, 'overlib.js'),
                    os.path.join(js_dir_path, 'overlib.js'))
    shutil.copyfile(
        os.path.join(js_path,
                     'otu_count_display.js'),
        os.path.join(js_dir_path,
                     'otu_count_display.js'))
    shutil.copyfile(os.path.join(js_path, 'jquery.js'),
                    os.path.join(js_dir_path, 'jquery.js'))
    shutil.copyfile(
        os.path.join(js_path,
                     'jquery.tablednd_0_5.js'),
        os.path.join(js_dir_path,
                     'jquery.tablednd_0_5.js'))

    # load tree for sorting OTUs
    ordered_otu_names = None
    if opts.tree is not None:
        try:
            f = open(opts.tree, 'U')
        except (TypeError, IOError):
            raise TreeMissingError(
                "Couldn't read tree file at path: %s" %
                tree_source)
        tree = parse_newick(f, PhyloNode)
        f.close()
        ordered_otu_names = [tip.Name for tip in tree.iterTips()]
    ordered_sample_names = None

    # load tree for sorting Samples
    if opts.sample_tree is not None:
        try:
            f = open(opts.sample_tree, 'U')
        except (TypeError, IOError):
            raise TreeMissingError(
                "Couldn't read tree file at path: %s" %
                tree_source)
        tree = parse_newick(f, PhyloNode)
        f.close()
        ordered_sample_names = [tip.Name for tip in tree.iterTips()]
    # if there's no sample tree, load sample map for sorting samples
    elif opts.map_fname is not None:
        lines = open(opts.map_fname, 'U').readlines()
        map = parse_mapping_file(lines)[0]
        ordered_sample_names = [row[0] for row in map]

    try:
        action = generate_heatmap_plots
    except NameError:
        action = None

    # Place this outside try/except so we don't mask NameError in action
    if action:
        action(
            num_otu_hits, otu_table, ordered_otu_names, ordered_sample_names,
            dir_path, js_dir_path, filename, fractional_values)

if __name__ == "__main__":
    main()
