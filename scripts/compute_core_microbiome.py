#!/usr/bin/env python
# File created on 12 Jun 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os.path import join

from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import xlim, ylim, xlabel, ylabel, plot, savefig
import numpy as np
from skbio.util import create_dir
from biom import load_table
from biom.exception import TableException

from qiime.util import (parse_command_line_parameters, make_option,
                        write_biom_table)
from qiime.core_microbiome import filter_table_to_core
from qiime.filter import sample_ids_from_metadata_description

script_info = {}
script_info['brief_description'] = "Identify the core microbiome."
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("",
     "Identify the core OTUs in otu_table.biom, defined as the OTUs that are present in at least 50% of the samples. Write the list of core OTUs to a text file, and a new BIOM file containing only the core OTUs.",
     "%prog -i otu_table.biom -o otu_table_core"))

script_info['script_usage'].append(
    ("",
     "Identify the core OTUs in otu_table.biom, defined as the OTUs that are present in all of the samples in the 'Fast' treatment (as specified in the mapping file). Write the list of core OTUs to a text file.",
     "%prog -i otu_table.biom -o otu_table_core_fast --mapping_fp map.txt --valid_states \"Treatment:Fast\""))


script_info['output_description'] = ""
script_info['required_options'] = [
    make_option(
        '-i',
        '--input_fp',
        type="existing_filepath",
        help='the input otu table in BIOM format'),
    make_option(
        '-o',
        '--output_dir',
        type="new_dirpath",
        help='directory to store output data'),
]
script_info['optional_options'] = [
    make_option('--max_fraction_for_core', type=float,
                help='the max fractions of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: %default]', default=1.0),
    make_option('--min_fraction_for_core', type=float,
                help='the min fractions of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: %default]', default=0.5),
    make_option('--num_fraction_for_core_steps', type=int,
                help='the number of evenly sizes steps to take between min_fraction_for_core and max_fraction_for_core [default: %default]', default=11),
    make_option('--otu_md', default='taxonomy',
                help='the otu metadata category to write to the output file [defualt: %default]'),
    make_option('--mapping_fp', type='existing_filepath',
                help='mapping file path (for use with --valid_states) [default: %default]'),
    make_option('--valid_states',
                help='description of sample ids to retain (for use with --mapping_fp) [default: %default]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    input_fp = opts.input_fp
    output_dir = opts.output_dir

    if opts.num_fraction_for_core_steps < 2:
        option_parser.error(
            "Must perform at least two steps. Increase --num_fraction_for_core_steps.")
    fractions_for_core = np.linspace(opts.min_fraction_for_core,
                                     opts.max_fraction_for_core,
                                     opts.num_fraction_for_core_steps)

    otu_md = opts.otu_md
    valid_states = opts.valid_states
    mapping_fp = opts.mapping_fp

    create_dir(output_dir)

    if valid_states and opts.mapping_fp:
        sample_ids = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'),
            valid_states)
        if len(sample_ids) < 1:
            option_parser.error(
                "--valid_states pattern didn't match any entries in mapping file: \"%s\"" %
                valid_states)
    else:
        # get core across all samples if user doesn't specify a subset of the
        # samples to work with
        sample_ids = None

    input_table = load_table(input_fp)

    otu_counts = []
    summary_figure_fp = join(output_dir, 'core_otu_size.pdf')
    for fraction_for_core in fractions_for_core:
        # build a string representation of the fraction as that gets used
        # several times
        fraction_for_core_str = "%1.0f" % (fraction_for_core * 100.)

        # prep output files
        output_fp = join(
            output_dir,
            'core_otus_%s.txt' %
            fraction_for_core_str)
        output_table_fp = join(
            output_dir,
            'core_table_%s.biom' %
            fraction_for_core_str)
        output_f = open(output_fp, 'w')

        try:
            core_table = filter_table_to_core(input_table,
                                              sample_ids,
                                              fraction_for_core)
        except TableException:
            output_f.write(
                "# No OTUs present in %s %% of samples." %
                fraction_for_core_str)
            output_f.close()
            otu_counts.append(0)
            continue

        # write some header information to file
        if sample_ids is None:
            output_f.write(
                "# Core OTUs across %s %% of samples.\n" %
                fraction_for_core_str)
        else:
            output_f.write(
                "# Core OTUs across %s %% of samples matching the sample metadata pattern \"%s\":\n# %s\n" %
                (fraction_for_core_str, valid_states, ' '.join(sample_ids)))

        # write the otu id and corresponding metadata for all core otus
        otu_count = 0
        for value, id_, md in core_table.iter(axis='observation'):
            output_f.write('%s\t%s\n' % (id_, md[otu_md]))
            otu_count += 1
        output_f.close()

        # write the core biom table
        write_biom_table(core_table, output_table_fp)

        # append the otu count to the list of counts
        otu_counts.append(otu_count)

    plot(fractions_for_core, otu_counts)
    xlim(min(fractions_for_core), max(fractions_for_core))
    ylim(0, max(otu_counts) + 1)
    xlabel(
        "Fraction of samples that OTU must be observed in to be considered 'core'")
    ylabel("Number of OTUs")
    savefig(summary_figure_fp)


if __name__ == "__main__":
    main()
