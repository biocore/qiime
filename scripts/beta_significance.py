#!/usr/bin/env python
# File created on 4 June 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Jose Antonio Navas", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.util import make_option
import os
from numpy import array

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
from cogent.maths.unifrac.fast_unifrac import fast_unifrac_permutations_file, TEST_ON_PAIRWISE, TEST_ON_TREE, TEST_ON_ENVS
from cogent.maths.unifrac.fast_unifrac import fast_p_test_file
from qiime.util import parse_command_line_parameters
from qiime.format import format_unifrac_sample_mapping
from biom import load_table


script_info = {}
script_info[
    'brief_description'] = "This script runs any of a set of common tests to determine if a sample is statistically significantly different from another sample"
script_info[
    'script_description'] = "The tests are conducted on each pair of samples present in the input otu table. See the unifrac tutorial online for more details (http://unifrac.colorado.edu/)"

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("Example:",
     "Perform 100 randomizations of sample/sequence assignments, and record the probability that sample 1 is phylogenetically different from sample 2, using the unifrac monte carlo significance test. The test is run for all pairs of samples.",
     "%prog -i otu_table.biom -t rep_set.tre -s unweighted_unifrac -o unw_sig.txt"))

script_info[
    'output_description'] = "The script outputs a tab delimited text file with each pair of samples and a p value representing the probability that a random sample/sequence assignment will result in more dissimilar samples than the actual pair of samples."

script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_filepath',
                help='input otu table in biom format'),
    make_option('-o', '--output_path', type='new_filepath',
                help='output results path'),
    make_option('-s', '--significance_test', type='choice',
                choices=[
                    'unweighted_unifrac',
                    'weighted_unifrac',
                    'weighted_normalized_unifrac',
                    'p-test'],
                help="significance test to use, options are 'unweighted_unifrac', 'weighted_unifrac', 'weighted_normalized_unifrac', or 'p-test'"),
    make_option(
        '-t',
        '--tree_path',
        type='existing_filepath',
        help='path to newick tree file'),

]
script_info['optional_options'] = [
    make_option('-n', '--num_iters', default=100, type="int",
                help='number of monte carlo randomizations [default: %default]'),
    make_option('-k', '--type_of_test', type='choice',
                choices=['all_together', 'each_pair', 'each_sample'], default='each_pair',
                help="type of significance test to perform, options are 'all_together', 'each_pair' or 'each_sample'. [default: %default]"),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    otu_table_fp = opts.input_path

    otu_table = load_table(otu_table_fp)

    sample_ids = otu_table.ids()
    otu_ids = otu_table.ids(axis='observation')

    # This is not memory safe: need to be able to load the otu table as ints
    otu_table_array = array(list(otu_table.iter_data(axis='observation')),
                            dtype='int')

    if opts.type_of_test == 'all_together':
        type_of_test = TEST_ON_TREE
        header_text = "sample\tp value\tp value (Bonferroni corrected)\n"
    elif opts.type_of_test == 'each_pair':
        type_of_test = TEST_ON_PAIRWISE
        header_text = "sample 1\tsample 2\tp value\tp value (Bonferroni corrected)\n"
    elif opts.type_of_test == 'each_sample':
        type_of_test = TEST_ON_ENVS
        header_text = "sample\tp value\tp value (Bonferroni corrected)\n"
        if opts.significance_test == 'p-test':
            raise RuntimeError(
                'significance test type "each_sample" not allowed for p-test')
    else:
        raise RuntimeError('significance test type "%s" not found' %
                           opts.type_of_test)

    # note, uses ugly temp file
    if opts.significance_test == 'unweighted_unifrac':
        tree_in = open(opts.tree_path, 'U')
        output_fp = opts.output_path + '_envs.tmp'

        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp, 'U')
        try:
            result = fast_unifrac_permutations_file(tree_in, envs_in,
                                                    weighted=False,
                                                    num_iters=opts.num_iters,
                                                    verbose=opts.verbose,
                                                    test_on=type_of_test)
        except ValueError as e:
            if e.message == ("No valid samples/environments found. Check"
                             " whether tree tips match otus/taxa present in"
                             " samples/environments"):
                raise ValueError(e.message + " and that the otu abundance is"
                                 " not relative.")
            raise e

        envs_in.close()
        os.remove(output_fp)

        of = open(opts.output_path, 'w')
        of.write("#unweighted unifrac significance test\n")
        of.write(header_text)
        for line in result:
            of.write('\t'.join(map(str, line)) + '\n')
        of.close()

    elif opts.significance_test == 'p-test':
        tree_in = open(opts.tree_path, 'U')
        output_fp = opts.output_path + '_envs.tmp'

        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp, 'U')

        result = fast_p_test_file(tree_in, envs_in,
                                  num_iters=opts.num_iters, verbose=opts.verbose, test_on=type_of_test)
        envs_in.close()
        os.remove(output_fp)
        of = open(opts.output_path, 'w')
        of.write(
            "#andy martin's p-test significance test\n")
        of.write(header_text)
        for line in result:
            of.write('\t'.join(map(str, line)) + '\n')
        of.close()

    elif opts.significance_test == 'weighted_unifrac':
        tree_in = open(opts.tree_path, 'U')
        output_fp = opts.output_path + '_envs.tmp'

        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp, 'U')

        result = fast_unifrac_permutations_file(tree_in, envs_in,
                                                weighted=True, num_iters=opts.num_iters, verbose=opts.verbose, test_on=type_of_test)
        envs_in.close()
        os.remove(output_fp)
        of = open(opts.output_path, 'w')
        of.write(
            "#weighted unifrac significance test\n")
        of.write(header_text)
        for line in result:
            of.write('\t'.join(map(str, line)) + '\n')
        of.close()

    elif opts.significance_test == 'weighted_normalized_unifrac':
        tree_in = open(opts.tree_path, 'U')
        output_fp = opts.output_path + '_envs.tmp'

        result = format_unifrac_sample_mapping(
            sample_ids, otu_ids, otu_table_array)
        of = open(output_fp, 'w')
        of.write('\n'.join(result))
        of.close()
        envs_in = open(output_fp, 'U')

        result = fast_unifrac_permutations_file(tree_in, envs_in,
                                                weighted='correct', num_iters=opts.num_iters, verbose=opts.verbose, test_on=type_of_test)
        envs_in.close()
        os.remove(output_fp)
        of = open(opts.output_path, 'w')
        of.write(
            "#weighted normalized unifrac significance test\n")
        of.write(
            "sample 1\tsample 2\tp value\tp value (Bonferroni corrected)\n")
        for line in result:
            of.write('\t'.join(map(str, line)) + '\n')
        of.close()

    else:
        raise RuntimeError('significance test "%s" not found' %
                           opts.significance_test)
if __name__ == "__main__":
    main()
