#!/usr/bin/env python
# File created on 09 Feb 2010
# file make_otu_heatmap.py

from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = [
    "Dan Knights",
    "Jose Carlos Clemente Litran",
    "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"


from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.util import MissingFileError
from qiime.make_otu_heatmap_html import get_otu_counts, \
    filter_by_otu_hits
import shutil
import os
from os.path import join
from qiime.util import get_qiime_project_dir
from qiime.parse import parse_mapping_file
from numpy import array, arange, asarray
#from qiime.parse import parse_otu_table
from qiime.make_otu_heatmap import plot_heatmap, \
    get_clusters, make_otu_labels, extract_metadata_column, \
    get_order_from_categories, get_order_from_tree, names_to_indices, \
    get_log_transform, get_overlapping_samples
from biom.parse import parse_biom_table

options_lookup = get_options_lookup()

# make_otu_heatmap_html.py
script_info = {}
script_info['brief_description'] = """Make heatmap of OTU table"""
script_info[
    'script_description'] = """Once an OTU table has been generated, it can be visualized using a heatmap. In these heatmaps each row corresponds to an OTU, and each column corresponds to a sample. The higher the relative abundance of an OTU in a sample, the more intense the color at the corresponsing position in the heatmap. By default, the OTUs (rows) will be clustered by UPGMA hierarchical clustering, and the samples (columns) will be presented in the order in which they appear in the OTU table. Alternatively, the user may pass in a tree to sort the OTUs (rows) or samples (columns), or both. For samples, the user may also pass in a mapping file. If the user passes in a mapping file and a metadata category, samples (columns in the heatmap) will be grouped by category value and subsequently clustered within each group."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Examples:""",
     """Using default values:""",
     """%prog -i otu_table.biom"""))
script_info['script_usage'].append(
    ("",
     """Different output directory (i.e., "otu_heatmap"):""",
     """%prog -i otu_table.biom -o otu_heatmap"""))
script_info['script_usage'].append(
    ("",
     """Sort the heatmap columns by the order in a mapping file, as follows:""",
     """%prog -i otu_table.biom -o otu_heatmap -m mapping_file.txt"""))
script_info['script_usage'].append(
    ("",
     """Sort the heatmap columns by Sample ID's and the heatmap rows by the order of tips in the tree, you can supply a tree as follows:""",
     """%prog -i otu_table.biom -o otu_heatmap -m mapping_file.txt -t tree_file.txt"""))
script_info['script_usage'].append(
    ("",
     """Group the heatmap columns by metadata category (e.g., GENDER), then cluster within each group:""",
     """%prog -i otu_table.biom -o otu_heatmap -m mapping_file.txt -c 'GENDER'"""))
script_info[
    'output_description'] = """The heatmap image is located in the specified output directory. It is formatted as a PDF file."""
script_info['required_options'] = [
    options_lookup['otu_table_as_primary_input']
]
script_info['optional_options'] = [
    options_lookup['output_dir'],
    make_option('-t', '--otu_tree', type='existing_filepath',
                help='Tree file to be used for sorting OTUs \
in the heatmap', default=None),
    make_option(
        '-m', '--map_fname', dest='map_fname', type='existing_filepath',
        help='Metadata mapping file to be used for sorting Samples in the \
heatmap.', default=None),
    make_option('-c', '--category', dest='category', type="string",
                help='Metadata category for sorting samples. Samples will be clustered within each category level using euclidean UPGMA.', default=None),
    make_option(
        '-s', '--sample_tree', dest='sample_tree', type='existing_filepath',
        help='Tree file to be used for sorting samples (e.g, output from \
upgma_cluster.py). If both this and the \
sample mapping file are provided, the mapping file is ignored.', default=None),
    make_option('--no_log_transform', action="store_true",
                help='Data will not be log-transformed. Without this option, \
all zeros will be set to a small \
value (default is 1/2 the smallest non-zero entry). Data will be translated \
to be non-negative after log transform, and num_otu_hits will be set to 0.',
                default=False),
    make_option('--suppress_row_clustering', action="store_true",
                help='No UPGMA clustering of OTUs (rows) is performed. If --otu_tree is provided, \
this flag is ignored.', default=False),
    make_option('--suppress_column_clustering', action="store_true",
                help='No UPGMA clustering of Samples (columns) is performed. If --map_fname is provided, \
this flag is ignored.', default=False),
    make_option('--absolute_abundance', action="store_true",
                help='Do not normalize samples to sum to 1.[default %default]', default=False),
    make_option('--log_eps', type="float",
                help='Small value to replace zeros for log transform. \
[default: 1/2 the smallest non-zero entry].', default=None),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    otu_table = parse_biom_table(open(opts.otu_table_fp, 'U'))
    lineages = []
    if (otu_table.ObservationMetadata is None or 'taxonomy' not in otu_table.ObservationMetadata[0]):
        print '\n\nWarning: The lineages are missing from the OTU table. If you used single_rarefaction.py to create your otu_table, make sure you included the OTU lineages.\n'
        lineages = [''] * len(otu_table.ObservationIds)
    else:
        for val, id, meta in otu_table.iterObservations():
            lineages.append([v for v in meta['taxonomy']])

    otu_labels = make_otu_labels(otu_table.ObservationIds, lineages)

    # Convert to relative abundance if requested
    if not opts.absolute_abundance:
        otu_table = otu_table.normObservationBySample()

    # Get log transform if requested
    if not opts.no_log_transform:
        if not opts.log_eps is None and opts.log_eps <= 0:
            print "Parameter 'log_eps' must be positive. Value was", opts.log_eps
            exit(1)
        otu_table = get_log_transform(otu_table, opts.log_eps)

    if opts.output_dir:
        if os.path.exists(opts.output_dir):
            dir_path = opts.output_dir
        else:
            try:
                os.mkdir(opts.output_dir)
                dir_path = opts.output_dir
            except OSError:
                pass
    else:
        dir_path = './'

    # Re-order samples by tree if provided
    if not opts.sample_tree is None:
        sample_order = get_order_from_tree(otu_table.SampleIds,
                                           open(opts.sample_tree, 'U'))

    # if there's no sample tree, sort samples by mapping file
    elif not opts.map_fname is None:
        lines = open(opts.map_fname, 'U').readlines()
        metadata = list(parse_mapping_file(lines))
        new_map, otu_table = get_overlapping_samples(metadata[0], otu_table)
        metadata[0] = new_map
        map_sample_ids = zip(*metadata[0])[0]

        # if there's a category, do clustering within each category
        if not opts.category is None:
            category_labels = \
                extract_metadata_column(otu_table.SampleIds,
                                        metadata, opts.category)
            sample_order = \
                get_order_from_categories(otu_table, category_labels)
        # else: just use the mapping file order
        else:
            ordered_sample_ids = []
            for sample_id in map_sample_ids:
                if sample_id in otu_table.SampleIds:
                    ordered_sample_ids.append(sample_id)
            sample_order = names_to_indices(
                otu_table.SampleIds,
                ordered_sample_ids)
    # if no tree or mapping file, perform upgma euclidean
    elif not opts.suppress_column_clustering:
        data = asarray([i for i in otu_table.iterObservationData()])
        sample_order = get_clusters(data, axis='column')
    # else just use OTU table ordering
    else:
        sample_order = arange(len(otu_table.SampleIds))

    # re-order OTUs by tree (if provided), or clustering
    if not opts.otu_tree is None:
        # open tree file
        try:
            f = open(opts.otu_tree, 'U')
        except (TypeError, IOError):
            raise MissingFileError(
                "Couldn't read tree file at path: %s" %
                opts.otu_tree)
        otu_order = get_order_from_tree(otu_table.ObservationIds, f)
        f.close()
    # if no tree or mapping file, perform upgma euclidean
    elif not opts.suppress_row_clustering:
        data = asarray([i for i in otu_table.iterObservationData()])
        otu_order = get_clusters(data, axis='row')
    # else just use OTU table ordering
    else:
        otu_order = arange(len(otu_table.ObservationIds))

    # otu_order and sample_order should be ids, rather than indices
    #  to use in sortObservationOrder/sortSampleOrder
    otu_id_order = [otu_table.ObservationIds[i] for i in otu_order]
    sample_id_order = [otu_table.SampleIds[i] for i in sample_order]

    # Re-order otu table, sampleids, etc. as necessary
    otu_table = otu_table.sortObservationOrder(otu_id_order)
    # otu_ids not used after: tagged for deletion
    otu_ids = array(otu_table.ObservationIds)[otu_order]
    otu_labels = array(otu_labels)[otu_order]
    otu_table = otu_table.sortSampleOrder(sample_id_order)
    sample_ids = array(otu_table.SampleIds)[sample_order]

    plot_heatmap(otu_table, otu_labels, sample_ids,
                 filename=join(dir_path, 'heatmap.pdf'))


if __name__ == "__main__":
    main()
