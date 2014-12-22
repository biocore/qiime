#!/usr/bin/env python
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2014, The QIIME project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from sys import exit
from tempfile import TemporaryFile

from numpy import nan
from biom import load_table

from qiime.index import compute_index
from qiime.util import parse_command_line_parameters, make_option, MetadataMap


script_info = {}
script_info[
    'brief_description'] = "Compute a taxonomy-based index"
script_info[
    'script_description'] = ("Compute a taxonomy-based index to score each "
                             "sample in a table. This method is based on the "
                             "microbial dysbiosis index constructed in Gevers "
                             "et al. 2014 "
                             "(http://www.ncbi.nlm.nih.gov/pubmed/24629344).")
script_info[
    'script_usage'] = [("Example:", "Compute the MD-index",
                        "%prog -i $PWD/table.biom.gz -e md -o table.md.txt"),
                       ("Example:", ("Compute the MD-index and add it to an "
                                     "existing mapping file"),
                        ("%prog -i $PWD/table.biom.gz -e md -o mapping.md.txt "
                         "-m map.txt")),
                       ("Example:", ("Compute an arbitrary index, where the "
                                     "abundance of Firmicutes and "
                                     "Fusobacteria is considered to be an "
                                     "increase for the index and the "
                                     "abundance of Bacteroidetes is "
                                     "considered to be a decrease. For "
                                     "instance, in Gevers et al 2014, the "
                                     "MD-index was formed by the organisms "
                                     "thought to be protective and "
                                     "inflammatory. The organisms increased "
                                     "abundance in healthy controls and "
                                     "protective were used as the increased "
                                     "abundance set, and the organisms "
                                     "thought to be inflammatory and at a "
                                     "decreased abundance in healthy controls "
                                     "were used as the decreased set."),
                        ("%prog -i $PWD/table.biom.gz --increased "
                         "p__Firmicutes,p__Fusobacteria --decreased "
                         "p__Bacteroidetes -o table.custom.txt"))]
script_info[
    'output_description'] = ("Two columns, the first being the sample ID and "
                             "the second being the index value for the "
                             "sample. If a mapping file is provided, then the "
                             "index values are merged into it.")
script_info['optional_options'] = [
    make_option('-i', '--input', type="existing_filepath",
                help='The input BIOM table'),
    make_option('-o', '--output', type='new_filepath',
                help='Path to where the output will be written'),
    make_option('--increased', type=str,
                help="Comma separated list of taxa considered to increase",
                default=None),
    make_option('--decreased', type=str,
                help="Comma separated list of taxa considered to decrease",
                default=None),
    make_option('-e', '--index', type=str,
                help="An existing index",
                default=None),
    make_option('-n', '--name', type=str,
                help="index name to use in the output",
                default=None),
    make_option('-m', '--mapping_file', type="existing_filepath",
                help="A mapping file to add the computed index to"),
    make_option('-k', '--key', type=str,
                help="Metadata key to use for computing [default: 'taxonomy']",
                default='taxonomy'),
    make_option('-s', '--show_indices', action='store_true',
                help="List known indices",
                default=False)]

script_info['version'] = __version__


known_indices = {
    'md': {'name': 'Microbial dysbiosis index',
           'source': 'http://www.ncbi.nlm.nih.gov/pubmed/24629344',
           'decreased': set(['g__Akkermansia',
                             'g__Faecalibacterium',
                             'g__Ruminococcus',
                             'g__Bacteroides',
                             'g__Parabacteroides']),
           'increased': set(['g__Escherichia',
                             'g__Fusobacterium',
                             'g__Haemophilus',
                             'g__Aggregatibacter',
                             'g__Sutterella',
                             'g__Veillonella'])}}


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    if opts.show_indices:
        for idx_key in sorted(known_indices):
            idx = known_indices[idx_key]
            print "%s: %s, %s" % (idx_key, idx['name'], idx['source'])
            print '\t', 'increased:'
            print '\n'.join(['\t\t%s' % t for t in idx['increased']])
            print '\t', 'decreased:'
            print '\n'.join(['\t\t%s' % t for t in idx['decreased']])
        exit(0)

    if opts.index is not None and known_indices.get(opts.index) is None:
        option_parser.error("%s is not a known index. Known indices are: %s"
                            % (opts.index, ','.join(known_indices.keys())))

    if opts.index is not None and (opts.increased or opts.decreased):
        option_parser.error("Cannot specify both an existing and custom index")

    if opts.index is None and opts.increased is None and \
            opts.decreased is None:
        option_parser.error("Must specify an existing or custom index")

    if opts.increased and opts.decreased is None:
        option_parser.error("Must specify decreased taxa")

    if opts.decreased and opts.increased is None:
        option_parser.error("Must specify increased taxa")

    if opts.index is not None:
        name = opts.name if opts.name is not None else opts.index
        increased = known_indices[opts.index]['increased']
        decreased = known_indices[opts.index]['decreased']
    else:
        name = opts.name if opts.name is not None else 'index'
        increased = set(opts.increased.split(','))
        decreased = set(opts.decreased.split(','))

    if opts.input is None:
        option_parser.error("Input not specified")

    if opts.output is None:
        option_parser.error("Output not specified")

    table = load_table(opts.input)

    if opts.mapping_file:
        mapping_file = open(opts.mapping_file, 'U')
        output_fp = TemporaryFile()
    else:
        mapping_file = None
        output_fp = open(opts.output, 'w')

    output_fp.write("#SampleID\t%s\n" % name)
    for id_, value in compute_index(table, increased, decreased, opts.key):
        output_fp.write("%s\t%f\n" % (id_, value))

    if opts.mapping_file:
        output_fp.seek(0)
        mapping_data = MetadataMap.mergeMappingFiles([output_fp, mapping_file],
                                                     no_data_value=nan)
        with open(opts.output, 'w') as fp:
            fp.write(str(mapping_data))

    output_fp.close()


if __name__ == "__main__":
    main()
