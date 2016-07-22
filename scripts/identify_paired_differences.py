#!/usr/bin/env python
# File created on 19 Jun 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Greg Caporaso", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from biom import load_table
from qiime.group import (
    extract_per_individual_state_metadata_from_sample_metadata,
    extract_per_individual_state_metadata_from_sample_metadata_and_biom)
from qiime.parse import parse_mapping_file_to_dict
from qiime.util import (parse_command_line_parameters,
                        make_option)
from qiime.filter import sample_ids_from_metadata_description
from qiime.stats import paired_difference_analyses

script_info = {}
script_info[
    'brief_description'] = "Generate plots and stats to test for change in some data point(s) with a state change on a per-individual basis."
script_info[
    'script_description'] = "This script provides a framework for paired-difference testing (i.e., analysis of data generated under a pre/post experimental design). In a pre/post experimental design, individuals are sampled before and after some 'treatment'. This code plots differences in values in the sample metadata (i.e., the mapping file) or observation counts in a BIOM table, and runs a (Bonferroni-corrected) one sample t-test on each sample metadata category or BIOM observation to determine if the mean of each distribution of pre/post differences differs from zero. If 'None' appears for the t score and p-values, this often means that the distribution of differences contained no variance, so the t-test could not be run. This can happen, for example, if the value passed for --valid_states is so restrictive that only a single sample is retained for analysis."
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Generate plots and stats for one category from the mapping file where the y-axis should be consistent across plots and the lines in the plots should be light blue.",
     "",
     "%prog -m map.txt --metadata_categories 'Streptococcus Abundance' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_results --ymin 0 --ymax 60 --line_color '#eeefff'"))
script_info['script_usage'].append(
    ("Generate plots and stats for three categories from the mapping file.",
     "",
     "%prog -m map.txt --metadata_categories 'Streptococcus Abundance,Phylogenetic Diversity,Observed OTUs' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_and_alpha_results"))
script_info['script_usage'].append(
    ("Generate plots for all observations in a biom file",
     "",
     "%prog -m map.txt -b otu_table.biom --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o otu_results"))
script_info['script_usage'].append(
    ("Generate plots for all observations in a biom file, but only including samples from individuals whose 'TreatmentResponse' was 'Improved' (as defined in the mapping file).",
     "",
     "%prog -m map.txt -b otu_table.biom --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o otu_results_improved_only --valid_states TreatmentResponse:Improved"))

script_info[
    'output_description'] = "The output of this script is plots of pre/post differences and associated statistics."

script_info['required_options'] = [
    make_option(
        '-m',
        '--mapping_fp',
        type="existing_filepath",
        help='the input metadata map filepath'),
    make_option(
        '-o',
        '--output_dir',
        type="new_filepath",
        help='directory where output files should be saved'),
    make_option(
        '-t',
        '--state_category',
        help='the mapping file column name to plot change over (usually has values like "pre-treatment" and "post-treatment")'),
    make_option(
        '-x',
        '--state_values',
        help='ordered list of state values to test change over (defines direction of graphs, generally something like "pre-treatment,post-treatment"). currently limited to two states.'),
    make_option(
        '-c',
        '--individual_id_category',
        help='the mapping file column name containing each individual\'s identifier (usually something like "personal_identifier")'),
]

script_info['optional_options'] = [
    make_option(
        '--ymin',
        default=None,
        type='float',
        help='set the minimum y-value across plots [default: determined on a per-plot basis]'),
    make_option(
        '--ymax',
        default=None,
        type='float',
        help='set the maximum y-value across plots [default: determined on a per-plot basis]'),
    make_option(
        '--metadata_categories',
        help='ordered list of the mapping file column names to test for paired differences (usually something like "StreptococcusAbundance,Phylogenetic Diversity") [default: %default]',
        default=None),
    make_option(
        '--observation_ids',
        help='ordered list of the observation ids to test for paired differences if a biom table is provided (usually something like "otu1,otu2") [default: compute paired differences for all observation ids]',
        default=None),
    make_option(
        '-b',
        '--biom_table_fp',
        help='path to biom table to use for computing paired differences [default: %default]',
        type='existing_filepath',
        default=None),
    make_option(
        '-s',
        '--valid_states',
        help="string describing samples that should be included based on their metadata (e.g. 'TreatmentResponse:Improved') [default: all samples are included in analysis]",
        default=None),
    make_option(
        '--line_color',
        help="color of lines in plots, useful if generating multiple plots in different runs of this script to overlay on top of one another. these can be specified as matplotlib color names, or as html hex strings [default: %default]",
        default="black"),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    mapping_fp = opts.mapping_fp
    state_values = opts.state_values.split(',')
    metadata_categories = opts.metadata_categories
    state_category = opts.state_category
    individual_id_category = opts.individual_id_category
    output_dir = opts.output_dir
    biom_table_fp = opts.biom_table_fp
    observation_ids = opts.observation_ids
    if not observation_ids is None:
        observation_ids = observation_ids.split(',')
    valid_states = opts.valid_states
    ymin = opts.ymin
    ymax = opts.ymax
    line_color = opts.line_color

    # validate the input - currently only supports either biom data
    # or mapping file data. if useful in the future it shouldn't be too
    # hard to allow the user to provide both.
    if metadata_categories and biom_table_fp:
        option_parser.error(
            "Can only pass --metadata_categories or --biom_table_fp, not both.")
    elif not (metadata_categories or biom_table_fp):
        option_parser.error(
            "Must pass either --metadata_categories or --biom_table_fp.")
    else:
        pass

    # parse the mapping file to a dict
    mapping_data = parse_mapping_file_to_dict(open(mapping_fp, 'U'))[0]

    # currently only support for pre/post (ie, two-state) tests
    if len(state_values) != 2:
        option_parser.error(
            "Exactly two state_values must be passed separated by a comma.")

    # filter mapping_data, if requested
    if valid_states:
        sample_ids_to_keep = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'), valid_states)
        for sid in mapping_data.keys():
            if sid not in sample_ids_to_keep:
                del mapping_data[sid]

    if biom_table_fp:
        biom_table = load_table(biom_table_fp)
        analysis_categories = observation_ids or biom_table.ids(axis='observation')
        personal_ids_to_state_values = \
            extract_per_individual_state_metadata_from_sample_metadata_and_biom(
                mapping_data,
                biom_table,
                state_category,
                state_values,
                individual_id_category,
                observation_ids=analysis_categories)
    else:
        analysis_categories = metadata_categories.split(',')
        personal_ids_to_state_values = \
            extract_per_individual_state_metadata_from_sample_metadata(
                mapping_data,
                state_category,
                state_values,
                individual_id_category,
                analysis_categories)

    paired_difference_analyses(personal_ids_to_state_values,
                               analysis_categories,
                               state_values,
                               output_dir,
                               line_color=line_color,
                               ymin=ymin,
                               ymax=ymax)


if __name__ == "__main__":
    main()
