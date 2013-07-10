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

from os.path import join
from math import ceil
from numpy import median, mean
import matplotlib
from matplotlib.pyplot import subplots
from pylab import savefig
from cogent.util.misc import create_dir
from cogent.maths.stats.test import t_one_sample
from qiime.parse import extract_per_individual_state_metadata_from_mapping_f
from qiime.util import (parse_command_line_parameters, 
                  make_option)

script_info = {}
script_info['brief_description'] = "Generate plots which illustrate the change in some data point(s) with a state change on a per-individual basis."
script_info['script_description'] = ""
script_info['script_usage'] = []
script_info['script_usage'].append(("Generate plots for two categories from the mapping file.","","%prog -m map.txt --metadata_categories 'Streptococcus Abundance,Veillonella Abundance' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_results"))
script_info['script_usage'].append(("Generate plots for four categories from the mapping file, where the y-axes should be set on a per-plot basis.","","%prog -m map.txt --metadata_categories 'Streptococcus Abundance,Veillonella Abundance,Phylogenetic Diversity,Observed OTUs' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_and_alpha_results --suppress_share_y_axis"))

script_info['output_description']= ""

script_info['required_options'] = [
 make_option('-m','--mapping_fp',type="existing_filepath",help='the input filepath'),
 make_option('-o','--output_dir',type="new_filepath",help='directory where output files should be saved'),
 make_option('--metadata_categories',help='ordered list of the mapping file column names to test for paired differences (usually something like "StreptococcusAbundance,Phylogenetic Diversity")'),
 make_option('--state_category',help='the mapping file column name to plot change over (usually has values like "pre-treatment" and "post-treatment")'),
 make_option('--state_values',help='ordered list of state values to test change over (defines direction of graphs, generally something like "pre-treatment,post-treatment"). currently limited to two states.'),
 make_option('--individual_id_category',help='the mapping file column name containing each individual\'s identifier (usually something like "personal_identifier")'),
]

script_info['optional_options'] = [
 make_option('--ymin',type="float",default=None,
             help='minimum value for y-axis [default: determined automatically]'),
 make_option('--ymax',type="float",default=None,
             help='maximum value for y-axis [default: determined automatically]'),
 make_option('--suppress_share_y_axis',default=False,action='store_true',
             help='set the scale of the y-axes will be set on a per-plot basis, rather than made consistent across all subplots [default: %default]')
 ]

script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    mapping_fp = opts.mapping_fp
    state_values = opts.state_values.split(',')
    metadata_categories = opts.metadata_categories.split(',')
    state_category = opts.state_category
    individual_id_category = opts.individual_id_category
    
    if len(state_values) != 2:
        option_parser.error("Exactly two state_values must be passed separated by a comma.")
    
    num_cols = 3
    num_subplots = len(metadata_categories)
    num_rows = int(ceil(num_subplots / num_cols))
    num_unused_subplots = (num_rows * num_cols) - num_subplots
    
    # create the subplot grid
    f, splts = subplots(num_rows,
                        num_cols,
                        sharex=False,
                        sharey=not opts.suppress_share_y_axis)
    for i in range(num_cols,num_cols - num_unused_subplots,-1):
        # blank out un-used subplot spaces
        try:
            splts[-1][i-1].axis('off')
        except TypeError:
            splts[i-1].axis('off')
    x_values = range(len(state_values))
    all_y_values = []
    
    create_dir(opts.output_dir)
    
    paired_difference_output_fp = join(opts.output_dir,'paired_difference_comparisons.txt')
    paired_difference_output_f = open(paired_difference_output_fp,'w')
    paired_difference_output_f.write("#Metadata category\tMean difference\tMedian difference\tt one sample\tt one sample parametric p-value\n")
    
    plot_output_fp = join(opts.output_dir,'plots.pdf')
    
    for category_number, metadata_category in enumerate(metadata_categories):
        personal_ids_to_responses = \
         extract_per_individual_state_metadata_from_mapping_f(
                                     open(mapping_fp,'U'),
                                     state_category,
                                     state_values,
                                     individual_id_category,
                                     metadata_category)
        
        differences = []
        for e in personal_ids_to_responses.values():
            if None in e:
                # no data for some of the entries, so skip this pid
                continue
            else:
                differences.append(e[1] - e[0])
        
        t_one_sample_results = t_one_sample(differences)
        paired_difference_output_f.write('\t'.join([metadata_category,
                                        str(mean(differences)),
                                        str(median(differences)),
                                        str(t_one_sample_results[0]),
                                        str(t_one_sample_results[1])]))
        paired_difference_output_f.write('\n')

        row_num = int(category_number/num_cols)
        col_num = int(category_number % num_cols)
        try:
            current_subplot = splts[row_num][col_num]
        except TypeError:
            # there is only one row, so the plot is 
            # access only by column number
            current_subplot = splts[col_num]
        current_x_values = []
        current_y_values = []
        
        for pid, data in personal_ids_to_responses.items():
            if None in data:
                # no data for some of the entries, so skip this pid
                continue
            all_y_values.extend(data)            
            current_subplot.plot(x_values,
                                 data,
                                 "black",
                                 linewidth=0.5)

        current_subplot.set_title(metadata_category,size=8)
        current_subplot.set_xticks(range(len(state_values)))
        current_subplot.set_xticklabels(state_values,size=6)
        
        
    #     for row_num in range(num_rows):
    #         splt = splts[row_num][col_num]
    #         slope = median_slopes[row_num]
    #         splt.text(0.40*(len(treatment_values)-1),
    #                   0.85*max(all_y_values),
    #                   "d=%1.2f" % slope,
    #                   size=8)
    #         min_y = opts.ymin or min(all_y_values)
    #         max_y = opts.ymax or max(all_y_values)
    #         splt.set_ylim(min_y,max_y)
    #         yticks = map(int,[min_y,mean([min_y,max_y]),max_y])
    #         splt.set_yticks(yticks)
    #         splt.set_xticks(range(len(treatment_values)))
    #                 
    #     splts[-1][col_num].set_xticklabels(treatment_values,size=6)
    #     splts[0][col_num].set_title(plot_category,size=10)
    # 
    # for i,response_value in enumerate(response_values):
    #     splts[i][0].set_ylabel(response_value,size=8)
    #     splts[i][0].set_yticklabels(yticks,size=6)
    
    savefig(plot_output_fp)
    paired_difference_output_f.close()
        

if __name__ == "__main__":
    main()