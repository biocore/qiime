#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski",
               "Julia Goodrich", "Daniel McDonald", "Antonio Gonzalez Pena",
               "Jesse Stombaugh", "Jose Carlos Clemente Litran",
               "Greg Caporaso", "Jai Ram Rideout", "Adam Robbins-Pianka"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"

from os.path import split, splitext, join
from sys import stdout, stderr

import numpy as np
from biom.table import Table
from biom import load_table

from qiime.util import (parse_command_line_parameters, make_option,
                        get_options_lookup, create_dir, write_biom_table)
from qiime.summarize_taxa import make_summary, add_summary_mapping
from qiime.parse import parse_mapping_file
from qiime.format import (
    write_add_taxa_summary_mapping, format_add_taxa_summary_mapping)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = (
    "Summarize taxa and store results in a new table or appended to an "
    "existing mapping file.")
script_info['script_description'] = (
    "The summarize_taxa.py script provides summary information of the "
    "representation of taxonomic groups within each sample. It takes an OTU "
    "table that contains taxonomic information as input. The taxonomic level "
    "for which the summary information is provided is designated with the -L "
    "option. The meaning of this level will depend on the format of the taxon "
    "strings that are returned from the taxonomy assignment step. The "
    "taxonomy strings that are most useful are those that standardize the "
    "taxonomic level with the depth in the taxonomic strings. For instance, "
    "the Greengenes database uses the following levels: Level 1 = Kingdom "
    "(e.g Bacteria), Level 2 = Phylum (e.g Actinobacteria), Level 3 = Class "
    "(e.g Actinobacteria), Level 4 = Order (e.g Actinomycetales), Level 5 = "
    "Family (e.g Streptomycetaceae), Level 6 = Genus (e.g Streptomyces), "
    "Level 7 = Species (e.g mirabilis). "
    "By default, the relative abundance of each taxonomic "
    "group will be reported, but the raw counts can be returned if -a is "
    "passed.\n\nBy default, taxa summary tables will be output in both "
    "classic (tab-separated) and BIOM formats. The BIOM-formatted taxa "
    "summary tables can be used as input to other QIIME scripts that accept "
    "BIOM files.")

script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Examples:",
     "Summarize taxa based at taxonomic levels 2, 3, 4, 5, and 6, and write "
     "resulting taxa tables to the directory './tax'",
     "%prog -i otu_table.biom -o ./tax"))

script_info['script_usage'].append(
    ("Examples:",
     "Summarize taxa based at taxonomic levels 2, 3, 4, 5, and 6, and write "
     "resulting mapping files to the directory './tax'",
     "%prog -i otu_table.biom -o tax_mapping/ -m Fasting_Map.txt"))

script_info['output_description'] = (
    "There are two possible output formats depending on whether or not a "
    "mapping file is provided with the -m option. If a mapping file is not "
    "provided, a table is returned where the taxonomic groups are each in a "
    "row and there is a column for each sample. If a mapping file is "
    "provided, the summary information will be appended to this file. "
    "Specifically, a new column will be made for each taxonomic group to "
    "which the relative abundances or raw counts will be added to the "
    "existing rows for each sample. The addition of the taxonomic information "
    "to the mapping file allows for taxonomic coloration of Principal "
    "coordinates plots in Emperor. As described in the Emperor "
    "documentation, principal coordinates plots can be dynamically colored based on "
    "any of the metadata columns in the mapping file. Dynamic coloration of "
    "the plots by the relative abundances of each taxonomic group can help to "
    "distinguish which taxonomic groups are driving the clustering patterns.")

script_info['required_options'] = [
    make_option('-i', '--otu_table_fp', dest='otu_table_fp',
                help='Input OTU table filepath [REQUIRED]',
                type='existing_filepath'),
]
script_info['optional_options'] = [
    make_option('-L', '--level', default='2,3,4,5,6', type='string',
                help='Taxonomic level to summarize by. [default: %default]'),
    make_option('-m', '--mapping',
                help='Input metadata mapping filepath. If supplied, then the '
                'taxon information will be added to this file. This option is '
                'useful for coloring PCoA plots by taxon abundance or to '
                'perform statistical tests of taxon/mapping associations.',
                type='existing_filepath'),
    make_option('--md_identifier', default='taxonomy', type='string',
                help='the relevant observation metadata key '
                '[default: %default]'),
    make_option('--md_as_string', default=False, action='store_true',
                help='metadata is included as string [default: metadata is '
                'included as list]'),
    make_option('-d', '--delimiter', action='store', type='string',
                dest='delimiter', default=';',
                help='Delimiter separating taxonomy levels. '
                '[default: %default]'),
    make_option('-a', '--absolute_abundance', action='store_true',
                dest='absolute_abundance', default=False,
                help='If present, the absolute abundance of the lineage in '
                'each sample is reported. By default, this script uses '
                'relative abundance [default: %default]'),
    make_option('-l', '--lower_percentage', type='float', default=None,
                help='If present, OTUs having higher absolute abundance are '
                'trimmed. To remove OTUs that make up more than 5% of the '
                'total dataset you would pass 0.05. [default: %default]'),
    make_option('-u', '--upper_percentage', type='float', default=None,
                help='If present, OTUs having lower absolute abundance are '
                'trimmed. To remove the OTUs that makes up less than 45% of '
                'the total dataset you would pass 0.45. [default: %default]'),
    make_option('-t', '--transposed_output', action='store_true',
                dest='transposed_output', default=False,
                help='If present, the output will be written transposed from '
                'the regular output. This is helpful in cases when you want '
                'to use Site Painter to visualize your data '
                '[default: %default]'),
    options_lookup['output_dir'],
    make_option('--suppress_classic_table_output', action='store_true',
                default=False, help='If present, the classic (TSV) format '
                'taxon table will not be created in the output directory. '
                'This option is ignored if -m/--mapping is present '
                '[default: %default]'),
    make_option('--suppress_biom_table_output', action='store_true',
                default=False, help='If present, the BIOM-formatted taxon '
                'table will not be created in the output directory. This '
                'option is ignored if -m/--mapping is present '
                '[default: %default]')
]
script_info['option_label'] = {'otu_table_fp': 'OTU table filepath',
                               'output_fp': 'Output filepath',
                               'mapping': 'QIIME-formatted mapping filepath',
                               'level': 'Summarize level',
                               'delimiter': 'Taxonomic delimiter',
                               'absolute_abundance': 'Use absolute abundance',
                               'lower_percentage': 'Top % of OTUs to remove',
                               'upper_percentage': 'Bottom % of OTUs to '
                                                   'remove'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    lower_percentage = opts.lower_percentage
    upper_percentage = opts.upper_percentage
    otu_table_fp = opts.otu_table_fp
    otu_table = load_table(otu_table_fp)
    delimiter = opts.delimiter
    mapping_fp = opts.mapping
    md_as_string = opts.md_as_string
    md_identifier = opts.md_identifier
    levels = opts.level.split(',')
    suppress_classic_table_output = opts.suppress_classic_table_output
    suppress_biom_table_output = opts.suppress_biom_table_output

    if upper_percentage is not None and lower_percentage is not None:
        raise ValueError(
            "upper_percentage and lower_percentage are mutually exclusive")

    if upper_percentage is not None and lower_percentage is not None and \
            mapping:
        raise ValueError("upper_percentage and lower_percentage can not be "
                         "using with mapping file")

    if upper_percentage is not None and \
            (upper_percentage < 0 or upper_percentage > 1.0):
        raise ValueError('max_otu_percentage should be between 0.0 and 1.0')

    if lower_percentage is not None and \
            (lower_percentage < 0 or lower_percentage > 1.0):
        raise ValueError('lower_percentage should be between 0.0 and 1.0')

    if mapping_fp:
        mapping_file = open(mapping_fp, 'U')
        mapping, header, comments = parse_mapping_file(mapping_file)

        # use the input Mapping file for producing the output filenames
        map_dir_path, map_fname = split(mapping_fp)
        map_basename, map_fname_ext = splitext(map_fname)
    else:
        if suppress_classic_table_output and suppress_biom_table_output:
            option_parser.error("Both classic and BIOM output formats were "
                                "suppressed.")

    if not opts.absolute_abundance:
        otu_table = otu_table.norm(axis='sample', inplace=False)

    # introduced output directory to will allow for multiple outputs
    if opts.output_dir:
        create_dir(opts.output_dir, False)
        output_dir_path = opts.output_dir
    else:
        output_dir_path = './'

    # use the input OTU table to produce the output filenames
    dir_path, fname = split(otu_table_fp)
    basename, fname_ext = splitext(fname)

    # Iterate over the levels and generate a summarized taxonomy for each
    for level in levels:
        if mapping_fp:
            # define output filename
            output_fname = join(output_dir_path,
                                map_basename + '_L%s.txt' % (level))

            summary, tax_order = add_summary_mapping(otu_table,
                                                     mapping,
                                                     int(level),
                                                     md_as_string,
                                                     md_identifier)

            write_add_taxa_summary_mapping(summary, tax_order, mapping,
                                           header, output_fname, delimiter)
        else:
            # define the output filename. The extension will be added to the
            # end depending on the output format
            output_fname = join(output_dir_path, basename + '_L%s' % level)

            summary, header = make_summary(otu_table,
                                           int(level),
                                           upper_percentage,
                                           lower_percentage,
                                           md_as_string,
                                           md_identifier)

            sample_ids = header[1:]

            observation_ids = []
            data = []
            for row in summary:
                # Join taxonomic levels to create an observation ID.
                observation_ids.append(delimiter.join(row[0]))
                data.append(row[1:])

            table = Table(np.asarray(data), observation_ids, sample_ids)
            if opts.transposed_output:
                table = table.transpose()

            if not suppress_classic_table_output:
                with open(output_fname + '.txt', 'w') as outfile:
                    outfile.write(table.to_tsv())

            if not suppress_biom_table_output:
                write_biom_table(table, output_fname + '.biom')


if __name__ == "__main__":
    main()
