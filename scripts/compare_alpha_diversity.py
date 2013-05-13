#!/usr/bin/env python
# File created on 06 Jun 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren", "Greg Caparaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, make_option, create_dir
from qiime.compare_alpha_diversity import (compare_alpha_diversities,
    _correct_compare_alpha_results, test_types, correction_types)
import os

script_info = {}
script_info['brief_description'] = """This script compares alpha diversities based on a two sample t-test using either parametric or non-parametric (Monte Carlo) methods."""
 
script_info['script_description'] = """
This script compares the alpha diversity of samples found in a collated alpha 
diversity file.  The comparison is done not between samples, but between groups
of samples. The groupings are created via the input 'category' passed with this 
script. Any samples which have the same value under the catgory will be grouped. 
For examples if your mapping file had a category called 'Treatment' that
separated your samples into three groups (Treatment='Control', Treatment='Drug',
Treatment='2xDose') passing 'Treatment' to this script would cause it to compare 
(Control,Drug), (Control,2xDose), (2xDose, Drug) alpha diversity values. 
By default the two sample t-test will be nonparametric (i.e. using Monte Carlo 
permutations to calculate the p-value), though the user has the option to make 
the test a parametric t-test.

The output is two tables located in the directory to be created: 
Comparisons.txt which has a format of comparison X (tval,pval) table where each 
row is a different group comparison, and the columns are a t-value or p-value 
(indicated by the header). Any iterations of a rarefaction at a given depth will 
be averaged. For instance, if your collated_alpha file had 10 iterations of the
rarefaction at depth 480, the scores for the alpha diversity metrics of those
10 iterations would be averaged (within sample). The iterations are not 
controlled by this script; when multiple_rarefactions.py is called, the -n option 
specifies the number of iterations that have occurred. The multiple comparison
correction takes into account the number of between group comparisons. If you do
not know the the rarefaction depth available or you want to use the deepest 
rarefaction level available then do not pass the -d or --depth option and it 
will default to using the deepest available. 
Averages.txt which has as format Treatmentgroup X average where the each row is
the average value of the alpha diversity of that sample group along with the
standard deviation in the averages. 

If t-scores and/or p-values are None for any of your comparisons there are three
possible reasons. The first is that there were undefined values in your collated
alpha diversity input file. This occurs if there were too few sequences in one 
or more of the samples in the groups involved in those comparisons to compute 
alpha diversity at that depth. You can either rerun %prog passing a lower 
value for --depth, or you can re-run alpha diversity after filtering samples
with too few sequences.
The second is that you had some comparison where each treatment was represented 
by only a single sample. It is not possible to perform a t_two_sample test on 
two samples each of length 1, so it will give you None,None for tval,pval 
instead.
The third possibility occurs when using the nonparamteric t test with small 
datasets where the monte carlo simulations return no p-value because the
distribution of the data has almost no variance.
The multiple comparisons correction will not penalize you for comparisons that
return as None regardless of origin. 

If the averages/standard deviations are None for any treatment group in the 
Averages.txt file the likely cause is that there is an \'n\\a\' value in the 
collated_alpha file that was passed. 
"""
 
script_info['script_usage'] = []

script_info['script_usage'].append(("Comparing alpha diversities",
"The following command takes the following input: a mapping file (which "
"associaties each sample with a number of characteristics), alpha diversity "
"metric (the results of collate_alpha for an alpha diverity metric, like "
"PD_whole_tree), depth (the rarefaction depth to use for comparison), "
"category (the category in the mapping file to determine which samples to "
"compare to each other), and output dir path (a path to the output directory to be created). A "
"nonparametric two sample t-test is run to compare the alpha diversities "
"using the default number of Monte Carlo permutations (999).",
"%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment -d 100 -o PD_d100.txt"))

script_info['script_usage'].append(("Parametric t-test",
"The following command runs a parametric two sample t-test using the "
"t-distribution instead of Monte Carlo permutations at rarefaction depth 100.",
"%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment -d 100 -o "
"PD_d100_parametric.txt -t parametric"))

script_info['script_usage'].append(("Parametric t-test",
"The following command runs a parametric two sample t-test using the "
"t-distribution instead of Monte Carlo permutations at the greatest depth available.",
"%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment -o "
"PD_dmax_parametric.txt -t parametric"))

script_info['output_description']= """
Script generates an output directory that is a table of TreatmentPair by (tval,pval). 
Each row corresponds to a comparison between two groups of treatment values. 
The columns are the tvals or pvals for that comparison.
"""

script_info['script_usage_output_to_remove'] = ['$PWD/PD_dmax_parametric.txt','$PWD/PD_d100_parametric.txt', '$PWD/PD_d100.txt']

script_info['required_options']=[
 make_option('-i',
  '--alpha_diversity_filepath',
  action='store',
  type='existing_filepath',
  dest='alpha_diversity_fp',
  help='path to collated alpha diversity file (as generated by '
       'collate_alpha.py) [REQUIRED]'),
 make_option('-m',
  '--mapping_filepath',
  action='store',
  type='existing_filepath',
  dest='mapping_fp',
  help='path to the mapping file [REQUIRED]'),
 make_option('-c',
  '--category',
  action='store',
  type='string',
  dest='category',
  help='category for comparison [REQUIRED]'),
 make_option('-o',
  '--output_fp',
  action='store',
  type='new_filepath',
  dest='output_fp',
  help='location of output file to be created [REQUIRED]')]

script_info['optional_options'] = [
 make_option('-t', '--test_type', type='choice', choices=test_types,
  help='the type of test to perform when calculating the p-values. Valid '
       'choices: ' + ', '.join(test_types) + '. If test_type is '
       'nonparametric, Monte Carlo permutations will be used to determine the '
       'p-value. If test_type is parametric, the num_permutations option will '
       'be ignored and the t-distribution will be used instead [default: '
       '%default]', default='nonparametric'),
 make_option('-n', '--num_permutations', type='int', default=999,
  help='the number of permutations to perform when calculating the '
       'p-value. Must be greater than 10. Only applies if test_type is '
       'nonparametric [default: %default]'),
  make_option('-p', '--correction_method', type='choice',
  choices=correction_types, help='method to use for correcting multiple '
  'comparisons. Available methods are bonferroni, fdr, or none. '
  '[default: %default]', default='bonferroni'),
   make_option('-d', '--depth', type='int', default=None, dest='depth',
  help='depth of rarefaction file to use [default: %default]')]


script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    if opts.num_permutations < 10:
        option_parser.error('Number of permuations must be greater than or '
                            'equal to 10.')

    rarefaction_lines = open(opts.alpha_diversity_fp, 'U')
    mapping_lines = open(opts.mapping_fp, 'U')
    category = opts.category
    depth = opts.depth

    ttest_result, alphadiv_avgs = compare_alpha_diversities(rarefaction_lines,
        mapping_lines, category, depth, opts.test_type, opts.num_permutations)
    
    rarefaction_lines.close()
    mapping_lines.close()

    corrected_result = _correct_compare_alpha_results(ttest_result,
        opts.correction_method)

    # write  results
    outfile = open(opts.output_fp, 'w')
    header = 'Comparison\tcat1_avg\tcat1_std\tcat2_avg\tcat2_std\tval\tpval'
    lines = [header]
    for k,v in corrected_result.items():
        t0, t1 = k.split(',') #if treatment values have ,s this will cause error
        lines.append('\t'.join(map(str,[k,alphadiv_avgs[t0][0],
            alphadiv_avgs[t0][1], alphadiv_avgs[t1][0],
            alphadiv_avgs[t1][1],v[0],v[1]])))
    outfile.write('\n'.join(lines))
    outfile.close()

if __name__ == "__main__":
    main()
