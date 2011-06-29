#!/usr/bin/env python
# File created on 4 June 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"


from qiime.util import make_option
import os

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
from cogent.maths.unifrac.fast_unifrac import fast_unifrac_permutations_file
from cogent.maths.unifrac.fast_unifrac import fast_p_test_file

from qiime.util import parse_command_line_parameters
from qiime.parse import parse_otu_table
from qiime.format import format_unifrac_sample_mapping


script_info = {}
script_info['brief_description'] = "This script runs any of a set of common tests to determine if a sample is statistically significantly different from another sample"
script_info['script_description'] = "The tests are conducted on each pair of samples present in the input otu table. See the unifrac tutorial online for more details (http://bmf2.colorado.edu/unifrac/tutorial.psp)"
script_info['script_usage'] = [("Example:","Perform 100 randomizations of sample/sequence assignments, and record the probability that sample 1 is phylogenetically different from sample 2, using the unifrac monte carlo significance test. The test is run for all pairs of samples.","%prog -i otu_table.txt -t rep_set.tre -s unweighted_unifrac -o unw_sig.txt")]
script_info['output_description']= "The script outputs a tab delimited text file with each pair of samples and a p value representing the probability that a random sample/sequence assignment will result in more dissimilar samples than the actual pair of samples."
script_info['required_options'] = [\
 make_option('-i', '--input_path',
     help='input otu table'), 
 make_option('-o', '--output_path',
     help='output results path'),
 make_option('-s', '--significance_test',
     help="significance test to use, options are 'unweighted_unifrac', 'weighted_unifrac', or 'p-test'"), 
 make_option('-t', '--tree_path',help='path to newick tree file'),

]
script_info['optional_options'] = [
 make_option('-n', '--num_iters', default=100, type="int",
     help='number of monte carlo randomizations'+\
     ' [default: %default]'),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    if opts.tree_path==None:
        if opts.significance_test in ['unweighted_unifrac','weighted_unifrac','p-test']:
            raise RuntimeError('please supply a phylogenetic tree for %s' %\
                opts.significance_test)

    # note, uses ugly temp file
    if opts.significance_test == 'unweighted_unifrac':
        tree_in = open(opts.tree_path,'U')
        otu_table_fp = opts.input_path
        output_fp = opts.output_path + '_envs.tmp'
    
        otu_table_lines = open(otu_table_fp, 'U')
        sample_ids, otu_ids, otu_table_array, lineages = \
        parse_otu_table(otu_table_lines)
        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp,'U')

        result = fast_unifrac_permutations_file(tree_in, envs_in,
            weighted=False, num_iters=opts.num_iters, verbose=opts.verbose)
        envs_in.close()
        os.remove(output_fp)
        
        of = open(opts.output_path,'w')
        of.write(\
            "#unweighted unifrac significance test\n")
        of.write(\
            "sample 1\tsample 2\tp value\tp value (Bonferroni corrected)\n")
        for line in result:
            of.write('\t'.join(map(str,line)) + '\n')
        of.close()

    elif opts.significance_test == 'p-test':
        tree_in = open(opts.tree_path,'U')
        otu_table_fp = opts.input_path
        output_fp = opts.output_path + '_envs.tmp'
    
        otu_table_lines = open(otu_table_fp, 'U')
        sample_ids, otu_ids, otu_table_array, lineages = \
        parse_otu_table(otu_table_lines)
        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp,'U')

        result = fast_p_test_file(tree_in, envs_in,
            num_iters=opts.num_iters, verbose=opts.verbose)
        envs_in.close()
        os.remove(output_fp)
        of = open(opts.output_path,'w')
        of.write(\
            "#andy martin's p-test significance test\n")
        of.write(\
            "sample 1\tsample 2\tp value\tp value (Bonferroni corrected)\n")
        for line in result:
            of.write('\t'.join(map(str,line)) + '\n')
        of.close()

    elif opts.significance_test == 'weighted_unifrac':
        tree_in = open(opts.tree_path,'U')
        otu_table_fp = opts.input_path
        output_fp = opts.output_path + '_envs.tmp'
    
        otu_table_lines = open(otu_table_fp, 'U')
        sample_ids, otu_ids, otu_table_array, lineages = \
        parse_otu_table(otu_table_lines)
        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp,'U')

        result = fast_unifrac_permutations_file(tree_in, envs_in,
            weighted=True, num_iters=opts.num_iters, verbose=opts.verbose)
        envs_in.close()
        os.remove(output_fp)
        of = open(opts.output_path,'w')
        of.write(\
            "#weighted unifrac significance test\n")
        of.write(\
            "sample 1\tsample 2\tp value\tp value (Bonferroni corrected)\n")
        for line in result:
            of.write('\t'.join(map(str,line)) + '\n')
        of.close()

    else:
        raise RuntimeError('significance test "%s" not found' %\
            opts.significance_test)
if __name__ == "__main__":
    main()
