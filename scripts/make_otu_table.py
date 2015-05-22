#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = [
    "Rob Knight",
    "Justin Kuczynski",
    "Jesse Stombaugh",
    "Greg Caporaso",
    "Adam Robbins-Pianka",
    "Sami Pietila"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from sys import argv, exit, stderr, stdout
from os.path import splitext

from qiime.filter import (get_seq_ids_from_seq_id_file,
                          get_seq_ids_from_fasta_file)
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, write_biom_table)
from qiime.parse import (parse_taxonomy, parse_mapping_file,
                         mapping_file_to_dict)
from qiime.make_otu_table import make_otu_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Make OTU table"""
script_info[
    'script_description'] = """The script make_otu_table.py tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU in the last column if a taxonomy file is supplied."""
script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Make OTU table",
     """Make an OTU table from an OTU map (i.e., result from pick_otus.py) and a taxonomy assignment file (i.e., result from assign_taxonomy.py). Write the output file to otu_table.biom.""",
     """%prog -i otu_map.txt -t tax_assignments.txt -o otu_table.biom"""))

script_info['script_usage'].append(
    ("Make OTU table, excluding OTU ids listed in a fasta file",
     """Make an OTU table, excluding the sequences listed in pynast_failures.fna. Note that the file pass as -e must end with either '.fasta' or '.fna'.""",
     "%prog -i otu_map.txt -o otu_table_no_pynast_failures.biom -e pynast_failures.fna"))

script_info['script_usage'].append(
    ("Make OTU table, excluding a list of OTU ids",
     """Make an OTU table, excluding the sequences listed in chimeric_seqs.txt""",
     "%prog -i otu_map.txt -o otu_table_non_chimeric.biom -e chimeric_seqs.txt"))

script_info['script_usage'].append(
    ("Make OTU table, passing a mapping file with sample metadata",
     """Make an OTU table from an OTU map (i.e., result from pick_otus.py). Write the output file to otu_table.biom.""",
     """%prog -i otu_map.txt -t tax_assignments.txt -o otu_table.biom -m mapping_file.txt"""))

script_info[
    'output_description'] = """The output of make_otu_table.py is a biom file, where the columns correspond to Samples and rows correspond to OTUs and the number of times a sample appears in a particular OTU."""

script_info['required_options'] = [
    options_lookup['otu_map_as_primary_input'],
    options_lookup['output_biom_fp'],
]

script_info['optional_options'] = [
    make_option(
        '-t', '--taxonomy', type='existing_filepath', dest='taxonomy_fname',
        help='Path to taxonomy assignment, containing the assignments of taxons to sequences (i.e., resulting txt file from assign_taxonomy.py) [default: %default]'),
    options_lookup['mapping_fp'],
    make_option('-e', '--exclude_otus_fp', type='existing_filepath',
                help=("path to a file listing OTU identifiers that should not be included in the "
                      "OTU table (e.g., the output of identify_chimeric_seqs.py) or a fasta "
                      "file where seq ids should be excluded (e.g., failures fasta file from align_seqs.py)"))
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    exclude_otus_fp = opts.exclude_otus_fp

    if not opts.taxonomy_fname:
        otu_to_taxonomy = None
    else:
        infile = open(opts.taxonomy_fname, 'U')
        otu_to_taxonomy = parse_taxonomy(infile)

    ids_to_exclude = []
    if exclude_otus_fp:
        if splitext(exclude_otus_fp)[1] in ('.fasta', '.fna'):
            ids_to_exclude = \
                get_seq_ids_from_fasta_file(open(exclude_otus_fp, 'U'))
        else:
            ids_to_exclude = \
                get_seq_ids_from_seq_id_file(open(exclude_otus_fp, 'U'))

    sample_metadata = None
    if opts.mapping_fp is not None:
        with open(opts.mapping_fp, 'U') as map_f:
            mapping_data, mapping_header, mapping_comments = \
                parse_mapping_file(map_f)

        sample_metadata = mapping_file_to_dict(mapping_data,
                                               mapping_header)

    with open(opts.otu_map_fp, 'U') as otu_map_f:
        biom_otu_table = make_otu_table(otu_map_f,
                                        otu_to_taxonomy=otu_to_taxonomy,
                                        otu_ids_to_exclude=ids_to_exclude,
                                        sample_metadata=sample_metadata)

    write_biom_table(biom_otu_table, opts.output_biom_fp)


if __name__ == "__main__":
    main()
