#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from os.path import join, split, splitext

from cogent.parse.tree import DndParser
from biom import load_table

from qiime.util import (add_filename_suffix, create_dir, get_options_lookup,
                        parse_command_line_parameters, make_option)
from qiime.simsam import simsam_range_to_files

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = ("Simulate samples for each sample in an OTU "
                            "table, using a phylogenetic tree.")
script_info['script_description'] = (
    "This script makes n samples related to each sample in an input otu table "
    "\n"
    "An input OTU table with 3 samples and n=2 will result in an output OTU "
    "table with 6 samples total: 3 clusters of 2 related samples."
    "\n"
    "To simulate each of the new samples, this script uses a sample in the "
    "input OTU table, and for each OTU in that sample the script traverses "
    "rootward on the tree a distance specified by '-d' to a point x. It then "
    "randomly selects a tip that decends from x, (call that new tip 'o2'), "
    "and reassigns all observations of the original OTU to the tip/OTU 'o2'.")
script_info['script_usage'] = [(
    "",
    "Create an OTU table with 3 related samples for each sample in "
    "otu_table.biom with dissimilarities of 0.001.",
    "%prog -i otu_table.biom -t rep_set.tre -o simsam_out1 -d .001 -n 3")]
script_info['script_usage'].append(
    ("",
     "Create OTU tables with 2, 3 and 4 related samples for each sample in "
     "otu_table.biom with dissimilarities of 0.001 and 0.01. Additionally "
     "create new mapping files with metadata for each of the new samples for "
     "use in downstream analyses.",
     "%prog -i otu_table.biom -t rep_set.tre -o simsam_out2 -d .001,.01 -n "
     "2,3,4 -m map.txt"))
script_info['output_description'] = """
The output directory will contain an OTU table with samples named:
'original_sample_0, original_sample_1 ...'

If a mapping file is provided via -m, an output mapping file containing the
replicated sample IDs (with all other metadata columns copied over) will also
be created.
"""
script_info['required_options'] = [
    make_option(
        '-i',
        '--otu_table',
        help='the input otu table',
        type='existing_filepath'),
    make_option(
        '-t',
        '--tree_file',
        help='tree file',
        type='existing_filepath'),
    options_lookup['output_dir'],
    make_option(
        '-d',
        '--dissim',
        help='dissimilarity between nodes up the tree, as a single value or '
        'comma-separated list of values'),
    make_option('-n', '--num', help='number of simulated samples per input '
                'sample, as a single value or comma-separated list of values')

]
script_info['optional_options'] = [
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help=('the mapping filepath. If provided, an output mapping '
                      'file containing the replicated sample IDs (with all '
                      'other metadata columns copied over) will also be '
                      'created [default: %default]'), default=None)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_dir = opts.output_dir
    create_dir(output_dir)

    otu_table_fp = opts.otu_table
    otu_table = load_table(otu_table_fp)

    tree_fh = open(opts.tree_file, 'U')
    tree = DndParser(tree_fh)
    tree_fh.close()

    mapping_fp = opts.mapping_fp
    if mapping_fp:
        mapping_f = open(mapping_fp, 'U')
        input_map_basename = splitext(split(mapping_fp)[1])[0]
    else:
        mapping_f = None
        input_map_basename = None

    input_table_basename = splitext(split(otu_table_fp)[1])[0]

    simsam_range_to_files(otu_table,
                          tree,
                          simulated_sample_sizes=map(int, opts.num.split(',')),
                          dissimilarities=map(float, opts.dissim.split(',')),
                          output_dir=output_dir,
                          mapping_f=mapping_f,
                          output_table_basename=input_table_basename,
                          output_map_basename=input_map_basename)


if __name__ == "__main__":
    main()
