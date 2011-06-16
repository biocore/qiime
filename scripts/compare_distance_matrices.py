#!/usr/bin/env python
# File created on 17 Mar 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from qiime.util import make_option
from cogent.maths.stats.test import mantel
from qiime.parse import parse_distmat
from qiime.format import format_p_value_for_num_iters
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,
                        make_compatible_distance_matrices)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Script for computing Mantel correlations between as set of distance matrices"
script_info['script_description'] = ""
script_info['script_usage'] = [("","Perform Mantel test on all pairs of four distance matrices, including 1000 Monte Carlo iterations. Write the output to mantel_out.txt.","mantel.py -i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt,weighted_unifrac_even100_dm.txt,unweighted_unifrac_even100_dm.txt -o mantel_out.txt -n 1000")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--input_dms',help='the input distance matrices, comma-separated'),\
 make_option('-o','--output_fp',help='the output filepath'),\
]
script_info['optional_options'] = [
 make_option('-n','--num_iterations',help='the number of iterations to perform',default=100,type='int')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    input_dm_fps = opts.input_dms.split(',')
    output_f = open(opts.output_fp,'w')
    output_f.write('DM1\tDM2\tNumber of entries\tMantel p-value\n')
    num_iterations = opts.num_iterations
    for i,fp1 in enumerate(input_dm_fps):
        for fp2 in input_dm_fps[i+1:]:
            (dm1_labels, dm1), (dm2_labels, dm2) =\
             make_compatible_distance_matrices(parse_distmat(open(fp1,'U')),
                                               parse_distmat(open(fp2,'U')))
            p = mantel(dm1,
                       dm2,
                       n=num_iterations)
            p_str = format_p_value_for_num_iters(p,num_iterations)
            output_f.write('%s\t%s\t%d\t%s\n' % (fp1,fp2,len(dm1_labels),p_str))
    output_f.close()

if __name__ == "__main__":
    main()