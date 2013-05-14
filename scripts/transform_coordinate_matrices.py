#!/usr/bin/env python
# File created on 21 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso","Justin Kuczynski","Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.util import make_option
from os.path import split, splitext, exists
from os import makedirs
from numpy import log10
from qiime.util import parse_command_line_parameters
from qiime.parse import fields_to_dict
from qiime.format import format_p_value_for_num_iters
from qiime.transform_coordinate_matrices import procrustes_monte_carlo,\
    get_procrustes_results

script_info={}
script_info['brief_description']="""Transform 2 coordinate matrices"""
script_info['script_description']="""This script transforms 2 coordinate matrices (e.g., the output of principal_coordinates.py) using procrustes analysis to minimize the distances between corresponding points. Monte Carlo simulations can additionally be performed (-r random trials are run) to estimate the probability of seeing an M^2 value as extreme as the actual M^2."""
script_info['script_usage']=[]
script_info['script_usage'].append(("Write the transformed procrustes matrices to file","","""%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt -o procrustes_output"""))

script_info['script_usage'].append(("Generate transformed procrustes matrices and monte carlo p-values","","""%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt -o mc_procrustes_output -r 1000""",))

script_info['output_description']="""Two transformed coordinate matrices corresponding to the two input coordinate matrices, and (if -r was specified) a text file summarizing the results of the Monte Carlo simulations."""
script_info['required_options']=[\
 make_option('-i','--input_fps',type='existing_filepaths',help='comma-separated input files'),\
 make_option('-o','--output_dir',type='new_dirpath',help='the output directory'),\
]
script_info['optional_options']=[\
 make_option('-r','--random_trials',type='int',
    help='Number of random permutations of matrix2 to perform. '+
    ' [default: (no Monte Carlo analysis performed)]',default=None),
 make_option('-d','--num_dimensions',type='int',default=3,
    help='Number of dimensions to include in output matrices'+
    ' [default: %default]'),
 make_option('-s','--sample_id_map_fp',
    type='existing_filepath',
    help='Map of original sample ids to new sample ids [default: %default]',
    default=None),
 make_option('--store_trial_details',
    help='Store PC matrices for individual trials [default: %default]',
    default=False,action='store_true'),
]

script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    random_trials = opts.random_trials
    if random_trials != None and random_trials < 10:
        option_parser.error('Must perform >= 10 trails for Monte Carlo analysis.')
        
    output_dir = opts.output_dir
    sample_id_map_fp = opts.sample_id_map_fp
    num_dimensions = opts.num_dimensions
    
    if not exists(output_dir): 
        makedirs(output_dir)
    
    if opts.store_trial_details:
        trial_output_dir = '%s/trial_details/' % output_dir
    else:
        trial_output_dir = None
  
    input_fp1 = opts.input_fps[0]
    input_fp2 = opts.input_fps[1]
    input_fp1_dir, input_fn1 = split(input_fp1)
    input_fp1_basename, input_fp1_ext = splitext(input_fn1)
    input_fp2_dir, input_fn2 = split(input_fp2)
    input_fp2_basename, input_fp2_ext = splitext(input_fn2)
    output_summary_fp = '%s/%s_%s_procrustes_results.txt' %\
     (output_dir,input_fp1_basename,input_fp2_basename)
    output_matrix1_fp = '%s/pc1_transformed.txt' % output_dir
    output_matrix2_fp = '%s/pc2_transformed.txt' % output_dir
    
    if sample_id_map_fp:
        sample_id_map = dict([(k,v[0]) \
         for k,v in fields_to_dict(open(sample_id_map_fp, "U")).items()])
    else:
        sample_id_map = None
    
    transformed_coords1, transformed_coords2, m_squared, randomized_coords2 =\
      get_procrustes_results(open(input_fp1,'U'),\
                             open(input_fp2,'U'),\
                             sample_id_map=sample_id_map,\
                             randomize=False,
                             max_dimensions=num_dimensions)
    output_matrix1_f = open(output_matrix1_fp,'w')
    output_matrix1_f.write(transformed_coords1)
    output_matrix1_f.close()
    output_matrix2_f = open(output_matrix2_fp,'w')
    output_matrix2_f.write(transformed_coords2)
    output_matrix2_f.close()
    
    if random_trials:
        summary_file_lines = ['FP1 FP2 Included_dimensions MC_p_value Count_better M^2']
        coords_f1 = list(open(input_fp1,'U'))
        coords_f2 = list(open(input_fp2,'U'))
        actual_m_squared, trial_m_squareds, count_better, mc_p_value =\
         procrustes_monte_carlo(coords_f1,\
                                coords_f2,\
                                trials=random_trials,\
                                max_dimensions=num_dimensions,
                                sample_id_map=sample_id_map,
                                trial_output_dir=trial_output_dir)
        # truncate the p-value to the correct number of significant
        # digits
        mc_p_value_str = format_p_value_for_num_iters(mc_p_value, random_trials)
        max_dims_str = str(num_dimensions or 'alldim')
        summary_file_lines.append('%s %s %s %s %d %1.3f' %\
         (input_fp1, input_fp2, str(max_dims_str), mc_p_value_str,\
          count_better, actual_m_squared))
        f = open(output_summary_fp,'w')
        f.write('\n'.join(summary_file_lines))
        f.write('\n')
        f.close()


if __name__ == "__main__":
    main()
