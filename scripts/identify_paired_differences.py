#!/usr/bin/env python
# File created on 19 Jun 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from biom.parse import parse_biom_table
from qiime.parse import (extract_per_individual_state_metadata_from_sample_metadata,
 extract_per_individual_state_metadata_from_sample_metadata_and_biom,
 parse_mapping_file_to_dict)
from qiime.util import (parse_command_line_parameters, 
                  make_option)
from qiime.filter import (filter_mapping_file_from_mapping_f,
                          sample_ids_from_metadata_description)
from qiime.stats import (run_paired_difference_analyses)

script_info = {}
script_info['brief_description'] = "Generate plots which illustrate the change in some data point(s) with a state change on a per-individual basis."
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("Generate plots for two categories from the mapping file.","","%prog -m map.txt --metadata_categories 'Streptococcus Abundance,Veillonella Abundance' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_results"))
script_info['script_usage'].append(("Generate plots for four categories from the mapping file, where the y-axes should be set on a per-plot basis.","","%prog -m map.txt --metadata_categories 'Streptococcus Abundance,Veillonella Abundance,Phylogenetic Diversity,Observed OTUs' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_and_alpha_results"))
script_info['script_usage'].append(("Generate plots for all observations in a biom file, where the y-axes should be set on a per-plot basis.","","%prog -m map.txt -b otu_table.biom --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o otu_results"))
script_info['script_usage'].append(("Generate plots for all observations in a biom file, where the y-axes should be set on a per-plot basis, and only including samples from individuals whose 'TreatmentResponse' was 'Improved' (as defined in the mapping file).","","%prog -m map.txt -b otu_table.biom --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o otu_results_improved_only --valid_states TreatmentResponse:Improved"))

script_info['output_description']= ""

script_info['required_options'] = [
 make_option('-m','--mapping_fp',type="existing_filepath",help='the input filepath'),
 make_option('-o','--output_dir',type="new_filepath",help='directory where output files should be saved'),
 make_option('--state_category',help='the mapping file column name to plot change over (usually has values like "pre-treatment" and "post-treatment")'),
 make_option('--state_values',help='ordered list of state values to test change over (defines direction of graphs, generally something like "pre-treatment,post-treatment"). currently limited to two states.'),
 make_option('--individual_id_category',help='the mapping file column name containing each individual\'s identifier (usually something like "personal_identifier")'),
]

script_info['optional_options'] = [
 make_option('--ymin',default=None,type='float',help='set the minimum y-value across plots [default: determined on a per-plot basis]'),
 make_option('--ymax',default=None,type='float',help='set the maximum y-value across plots [default: determined on a per-plot basis]'),
  make_option('--metadata_categories',help='ordered list of the mapping file column names to test for paired differences (usually something like "StreptococcusAbundance,Phylogenetic Diversity") [default: %default]',default=None),
  make_option('--observation_ids',help='ordered list of the observation ids to test for paired differences if a biom table is provided (usually something like "otu1,otu2") [default: compute paired differences for all observation ids]',default=None),
  make_option('-b','--biom_table_fp',help='path to biom table to use for computing paired differences [default: %default]', default=None),
  make_option('-s','--valid_states', help="string describing valid samples that should be included based on their metadata (e.g. 'TreatmentResponse:Improved') [default: %default]",default=None),
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
    
    # validate the input - currently only supports either biom data
    # or mapping file data. if useful in the future it shouldn't be too
    # hard to allow the user to provide both.
    if metadata_categories and biom_table_fp:
        option_parser.error("Can only pass --metadata_categories or --biom_table_fp, not both.")
    elif not (metadata_categories or biom_table_fp):
        option_parser.error("Must pass either --metadata_categories or --biom_table_fp.")
    else:
        pass
    
    # parse the mapping file to a dict
    mapping_data = parse_mapping_file_to_dict(open(mapping_fp,'U'))[0]
    
    # currently only support for pre/post tests
    if len(state_values) != 2:
        option_parser.error("Exactly two state_values must be passed separated by a comma.")
    
    # filter the mapping file, if requested
    if valid_states:
        sample_ids_to_keep = sample_ids_from_metadata_description(
                              open(mapping_fp,'U'),valid_states)
        mapping_f = \
         filter_mapping_file_from_mapping_f(
          open(mapping_fp,'U'), sample_ids_to_keep).split('\n')
        mapping_data = parse_mapping_file_to_dict(mapping_f)[0]
    
    if biom_table_fp:
        biom_table = parse_biom_table(open(biom_table_fp,'U'))
        analysis_categories = observation_ids or biom_table.ObservationIds
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

    run_paired_difference_analyses(personal_ids_to_state_values,
                                   analysis_categories,
                                   state_values,
                                   output_dir)

        

if __name__ == "__main__":
    main()