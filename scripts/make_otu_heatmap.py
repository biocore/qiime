#!/usr/bin/env python
from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = [
    "Dan Knights",
    "Jose Carlos Clemente Litran",
    "Yoshiki Vazquez Baeza",
    "Greg Caporaso",
    "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"


import numpy as np
from biom import load_table

from qiime.make_otu_heatmap import (
    plot_heatmap, get_clusters, make_otu_labels, extract_metadata_column,
    get_order_from_categories, get_order_from_tree, names_to_indices,
    get_log_transform, get_overlapping_samples)
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, MissingFileError)
from qiime.parse import parse_mapping_file

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Plot heatmap of OTU table"""
script_info['script_description'] = (
    "This script visualizes an OTU table as a heatmap where each row "
    "corresponds to an OTU and each column corresponds to a sample. The "
    "higher the relative abundance of an OTU in a sample, the more intense "
    "the color at the corresponsing position in the heatmap. By default, the "
    "OTUs (rows) will be clustered by UPGMA hierarchical clustering, and the "
    "samples (columns) will be presented in the order in which they appear in "
    "the OTU table. Alternatively, the user may supply a tree to sort the "
    "OTUs (rows) or samples (columns), or both. The user may also pass in a "
    "mapping file for sorting samples. If the user passes in a mapping file "
    "and a metadata category, samples (columns) will be grouped by category "
    "value and subsequently clustered within each group.")
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("",
     "Generate a heatmap as a PDF using all default values:",
     "%prog -i otu_table.biom -o heatmap.pdf"))
script_info['script_usage'].append(
    ("",
     "Generate a heatmap as a PNG:",
     "%prog -i otu_table.biom -o heatmap.png -g png"))
script_info['script_usage'].append(
    ("",
     "Sort the heatmap columns (samples) by the order of samples in the "
     "mapping file",
     "%prog -i otu_table.biom -o heatmap_sorted_samples.pdf -m "
     "mapping_file.txt"))
script_info['script_usage'].append(
    ("",
     "Sort the heatmap columns (samples) by the order of samples in the "
     "mapping file, and sort the heatmap rows by the order of tips in the "
     "tree:",
     "%prog -i otu_table.biom -o heatmap_sorted.pdf -m mapping_file.txt -t "
     "rep_set.tre"))
script_info['script_usage'].append(
    ("",
     "Group the heatmap columns (samples) by metadata category (e.g., "
     "Treatment), then cluster within each group:""",
     "%prog -i otu_table.biom -o heatmap_grouped_by_Treatment.pdf -m "
     "mapping_file.txt -c Treatment"))

script_info['output_description'] = (
    "A single output file is created containing the heatmap of the OTU "
    "table (a PDF file by default).")

script_info['required_options'] = [
    options_lookup['otu_table_as_primary_input'],
    options_lookup['output_fp']
]

script_info['optional_options'] = [
    make_option('-t', '--otu_tree', type='existing_filepath', help='Tree file '
                'to be used for sorting OTUs in the heatmap', default=None),
    make_option('-m', '--map_fname', dest='map_fname',
                type='existing_filepath', help='Metadata mapping file to be '
                'used for sorting Samples in the heatmap.', default=None),
    make_option('-c', '--category', dest='category', type="string",
                help='Metadata category for sorting samples. Samples will be '
                'clustered within each category level using euclidean UPGMA.',
                default=None),
    make_option('-s', '--sample_tree', dest='sample_tree',
                type='existing_filepath', help='Tree file to be used for '
                'sorting samples (e.g, output from upgma_cluster.py). If both '
                'this and the sample mapping file are provided, the mapping '
                'file is ignored.', default=None),
    make_option('-g', '--imagetype',
                help='type of image to produce (i.e. png, svg, pdf) '
                '[default: %default]', default='pdf', type="choice",
                choices=['pdf', 'png', 'svg']),
    make_option('--no_log_transform', action="store_true",
                help='Data will not be log-transformed. Without this option, '
                'all zeros will be set to a small value (default is 1/2 the '
                'smallest non-zero entry). Data will be translated to be '
                'non-negative after log transform, and num_otu_hits will be '
                'set to 0.', default=False),
    make_option('--suppress_row_clustering', action="store_true",
                help='No UPGMA clustering of OTUs (rows) is performed. If '
                '--otu_tree is provided, this flag is ignored.',
                default=False),
    make_option('--suppress_column_clustering', action="store_true",
                help='No UPGMA clustering of Samples (columns) is performed. '
                'If --map_fname is provided, this flag is ignored.',
                default=False),
    make_option('--absolute_abundance', action="store_true",
                help='Do not normalize samples to sum to 1 [default: %default]',
                default=False),
    make_option('--color_scheme', default="YlGn",
                help="color scheme for figure. see "
                     "http://matplotlib.org/examples/color/"
                     "colormaps_reference.html for choices "
                     "[default: %default]"),
    make_option('--width',
                help='width of the figure in inches [default: %default]',
                default=5, type='float'),
    make_option('--height',
                help='height of the figure in inches [default: %default]',
                default=5, type='float'),
    make_option('--dpi',
                help='resolution of the figure in dots per inch '
                     '[default: value of savefig.dpi in matplotlibrc file]',
                type='int', default=None),
    make_option('--obs_md_category', default="taxonomy",
                help="observation metadata category to plot "
                     "[default: %default]"),
    make_option('--obs_md_level', default=None, type="int",
                help="the level of observation metadata to plot for "
                     "hierarchical metadata [default: lowest level]")
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    otu_table = load_table(opts.otu_table_fp)
    obs_md_category = opts.obs_md_category
    obs_md_level = opts.obs_md_level
    if obs_md_level is None:
        # grab the last level if the user didn't specify a level
        obs_md_level = -1
    else:
        # convert to 0-based indexing
        obs_md_level -= 1
    obs_md = otu_table.metadata(axis='observation')

    obs_md_labels = []
    if (obs_md is None or obs_md_category not in obs_md[0]):
        obs_md_labels = [['']] * len(otu_table.ids(axis='observation'))
    else:
        for _, _, md in otu_table.iter(axis='observation'):
            current_md = md[obs_md_category]
            if obs_md_level < len(current_md):
                current_md_at_level = current_md[obs_md_level]
            else:
                current_md_at_level = ''
            obs_md_labels.append([current_md_at_level])

    otu_labels = make_otu_labels(otu_table.ids(axis='observation'),
                                 obs_md_labels)

    # Convert to relative abundance if requested
    if not opts.absolute_abundance:
        otu_table = otu_table.norm(axis='observation')

    # Get log transform if requested
    if not opts.no_log_transform:
        otu_table = get_log_transform(otu_table)

    # Re-order samples by tree if provided
    if opts.sample_tree is not None:
        sample_order = get_order_from_tree(otu_table.ids(),
                                           open(opts.sample_tree, 'U'))

    # if there's no sample tree, sort samples by mapping file
    elif opts.map_fname is not None:
        lines = open(opts.map_fname, 'U').readlines()
        metadata = list(parse_mapping_file(lines))
        new_map, otu_table = get_overlapping_samples(metadata[0], otu_table)
        metadata[0] = new_map
        map_sample_ids = zip(*metadata[0])[0]

        # if there's a category, do clustering within each category
        if opts.category is not None:
            category_labels = extract_metadata_column(otu_table.ids(),
                                                      metadata, opts.category)
            sample_order = get_order_from_categories(otu_table,
                                                     category_labels)
        # else: just use the mapping file order
        else:
            ordered_sample_ids = []
            for sample_id in map_sample_ids:
                if otu_table.exists(sample_id):
                    ordered_sample_ids.append(sample_id)
            sample_order = names_to_indices(
                otu_table.ids(),
                ordered_sample_ids)
    # if no tree or mapping file, perform upgma euclidean
    elif not opts.suppress_column_clustering:
        data = np.asarray([i for i in otu_table.iter_data(axis='observation')])
        sample_order = get_clusters(data, axis='column')
    # else just use OTU table ordering
    else:
        sample_order = np.arange(len(otu_table.ids()))

    # re-order OTUs by tree (if provided), or clustering
    if opts.otu_tree is not None:
        # open tree file
        try:
            f = open(opts.otu_tree, 'U')
        except (TypeError, IOError):
            raise MissingFileError("Couldn't read tree file at path: %s" %
                                   opts.otu_tree)
        otu_order = get_order_from_tree(otu_table.ids(axis='observation'), f)
        f.close()
    # if no tree or mapping file, perform upgma euclidean
    elif not opts.suppress_row_clustering:
        data = np.asarray([i for i in otu_table.iter_data(axis='observation')])
        otu_order = get_clusters(data, axis='row')
    # else just use OTU table ordering
    else:
        otu_order = np.arange(len(otu_table.ids(axis='observation')))

    # otu_order and sample_order should be ids, rather than indices
    #  to use in sortObservationOrder/sortSampleOrder
    otu_id_order = [otu_table.ids(axis='observation')[i] for i in otu_order]
    sample_id_order = [otu_table.ids()[i] for i in sample_order]

    # Re-order otu table, sampleids, etc. as necessary
    otu_table = otu_table.sort_order(otu_id_order, axis='observation')
    otu_labels = np.array(otu_labels)[otu_order]
    otu_table = otu_table.sort_order(sample_id_order)
    sample_labels = otu_table.ids()

    plot_heatmap(otu_table, otu_labels, sample_labels, opts.output_fp,
                 imagetype=opts.imagetype, width=opts.width,
                 height=opts.height, dpi=opts.dpi,
                 color_scheme=opts.color_scheme)


if __name__ == "__main__":
    main()
