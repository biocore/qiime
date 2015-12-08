#!/usr/bin/env python
from __future__ import division

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2014, The QIIME project"
__credits__ = ["Daniel McDonald", "Greg Caporaso", "Jai Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"

from sys import exit
from tempfile import TemporaryFile

from numpy import nan
from biom import load_table

from qiime.compute_taxonomy_ratios import compute_index
from qiime.util import parse_command_line_parameters, make_option, MetadataMap

known_indices = {
    'md': {'name': 'Microbial dysbiosis index',
    'source': 'http://www.ncbi.nlm.nih.gov/pubmed/24629344',
    'decreased': set(['g__Dialister',
    'g__Faecalibacterium',
    'g__Ruminococcus',
    'g__Sutterella',
    'g__Rikenellaceae',
    'g__Parabacteroides',
    'g__Bacteroides',
    'g__Lachnospiraceae',
    'g__Coprococcus',
    'g__Erysipelotrichaceae',
    'g__Dorea',
    'g__Ruminococcaceae',
    'g__Oscillospira',
    'g__Bilophila']),
    'increased': set(['g__Escherichia',
    'g__Haemophilus',
    'g__Fusobacterium',
    'g__Veillonella'])}}

script_info = {}
script_info[
    'brief_description'] = "Compute the log ratio abundance of specified taxonomic groups."
script_info[
    'script_description'] = script_info['brief_description'] + (" This method is based on the "
                             "microbial dysbiosis index described in Gevers "
                             "et al. 2014: "
                             "http://www.ncbi.nlm.nih.gov/pubmed/24629344")
script_info[
    'script_usage'] = [("Example:", "Compute the microbial dysbiosis (MD) index",
                        "%prog -i table.biom.gz -e md -o md.txt"),
                       ("Example:", ("Compute the microbial dysbiosis (MD) index and add it to an "
                                     "existing mapping file"),
                        ("%prog -i table.biom.gz -e md -o map_w_md.txt -m map.txt")),
                       ("Example:", (
                        "Compute the log of the abundance of p__Firmicutes "
                        "plus p__Fusobacteria divided by the abundance of "
                        "p__Bacteroidetes and write the results to "
                        "custom_index.txt."),
                        ("%prog -i table.biom.gz --increased "
                         "p__Firmicutes,p__Fusobacteria --decreased "
                         "p__Bacteroidetes -o custom_index.txt"))]
script_info[
    'output_description'] = ("By default, a minimal QIIME mapping file is "
                             "created containing two columns: SampleID and "
                             "the index. If -m is provided, the information "
                             "in that mapping file is merged into the default "
                             "output mapping file.")

script_info['required_options'] = []

script_info['optional_options'] = [
    make_option('-i', '--input', type="existing_filepath",
                help='The input BIOM table [REQUIRED if not passing -s]'),
    make_option('-o', '--output', type='new_filepath',
                help='Path to where the output will be written; this will be '
                     'a new sample metadata mapping file '
                     '[REQUIRED if not passing -s]'),
    make_option('--increased', type=str,
                help="Comma-separated list of taxa whose abundances are "
                     "included in the numerator of the ratio "
                     "[REQUIRED if not passing -s or -e]",
                default=None),
    make_option('--decreased', type=str,
                help="Comma-separated list of taxa whose abundances are "
                     "included in the denominator of the ratio "
                     "[REQUIRED if not passing -s or -e]",
                default=None),
    make_option('-e', '--index', type='choice', choices=known_indices.keys(),
                help=("Apply an existing index. Options are: %s "
                       "[REQUIRED if not passing -s or --increased and --decreased]"
                       % ','.join(known_indices.keys())),
                default=None),
    make_option('-n', '--name', type=str,
                help="Column name for the index in the output file [default: "
                     "'index', or value passed as -e if provided]",
                default=None),
    make_option('-m', '--mapping_file', type="existing_filepath",
                help="A mapping file containing data that should be included in the output file "
                     "[default: no additional mapping file data is included in output]"),
    make_option('-k', '--key', type=str,
                help="Metadata key to use for computing index [default: %default]",
                default='taxonomy'),
    make_option('-s', '--show_indices', action='store_true',
                help="List known indices and exit [default: %default]",
                default=False)]

script_info['version'] = __version__

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
        output_file = TemporaryFile()
    else:
        mapping_file = None
        output_file = open(opts.output, 'w')

    output_file.write("#SampleID\t%s\n" % name)
    for id_, value in compute_index(table, increased, decreased, opts.key):
        output_file.write("%s\t%f\n" % (id_, value))

    if opts.mapping_file:
        output_file.seek(0)
        mapping_data = MetadataMap.mergeMappingFiles([output_file, mapping_file],
                                                     no_data_value=nan)
        with open(opts.output, 'w') as f:
            f.write(str(mapping_data))

    output_file.close()


if __name__ == "__main__":
    main()
