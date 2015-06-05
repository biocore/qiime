#!/usr/bin/env python
from __future__ import division

__author__ = "Logan Knecht"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Logan Knecht", "Michael Dwan", "Damien Coy", "Jai Ram Rideout",
               "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import exists
from skbio.util import create_dir
from qiime.compare_categories import compare_categories, methods
from qiime.util import (parse_command_line_parameters, make_option,
                        get_options_lookup)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Analyzes statistical significance of sample groupings using distance matrices"""
script_info['script_description'] = """
This script allows for the analysis of the strength and statistical
significance of sample groupings using a distance matrix as the primary input.
Several statistical methods are available: adonis, ANOSIM, BIO-ENV, Moran's I,
MRPP, PERMANOVA, PERMDISP, and db-RDA.

Note: R's vegan and ape packages are used to compute many of these methods, and
for the ones that are not, their implementations are based on the
implementations found in those packages. It is recommended to read through the
detailed descriptions provided by the authors (they are not reproduced here)
and to refer to the primary literature for complete details, including the
methods' assumptions. To view the documentation of a method in R, prepend a
question mark before the method name. For example:

?vegan::adonis

The following are brief descriptions of the available methods:

adonis - Partitions a distance matrix among sources of variation in order to
describe the strength and significance that a categorical or continuous
variable has in determining variation of distances. This is a nonparametric
method and is nearly equivalent to db-RDA (see below) except when distance
matrices constructed with semi-metric or non-metric dissimilarities are
provided, which may result in negative eigenvalues. adonis is very similar to
PERMANOVA, though it is more robust in that it can accept either categorical or
continuous variables in the metadata mapping file, while PERMANOVA can only
accept categorical variables. See vegan::adonis for more details.

ANOSIM - Tests whether two or more groups of samples are significantly
different based on a categorical variable found in the metadata mapping file.
You can specify a category in the metadata mapping file to separate
samples into groups and then test whether there are significant differences
between those groups. For example, you might test whether 'Control' samples are
significantly different from 'Fast' samples. Since ANOSIM is nonparametric,
significance is determined through permutations. See vegan::anosim for more
details.

BIO-ENV - Finds subsets of variables whose Euclidean distances (after scaling
the variables) are maximally rank-correlated with the distance matrix. For
example, the distance matrix might contain UniFrac distances between
communities, and the variables might be numeric environmental variables (e.g.,
pH and latitude). Correlation between the community distance matrix and
Euclidean environmental distance matrix is computed using Spearman's rank
correlation coefficient (rho). This method will only accept categories that are
numerical (continuous or discrete). This is currently the only method in the
script that accepts more than one category (via -c). See vegan::bioenv for more
details. This method is also known as BEST (previously called BIO-ENV) in the
PRIMER-E software package.

Moran's I - This method uses the numerical (e.g. geographical) data supplied to
identify what type of spatial configuration occurs in the samples. For example,
are they dispersed, clustered, or of no distinctly noticeable configuration
when compared to each other? This method will only accept a category that is
numerical. See ape::Moran.I for more details.

MRPP - This method tests whether two or more groups of samples are
significantly different based on a categorical variable found in the metadata
mapping file. You can specify a category in the metadata mapping file to
separate samples into groups and then test whether there are significant
differences between those groups. For example, you might test whether 'Control'
samples are significantly different from 'Fast' samples. Since MRPP is
nonparametric, significance is determined through permutations. See
vegan::mrpp for more details.

PERMANOVA - This method is very similar to adonis except that it only accepts a
categorical variable in the metadata mapping file. It uses an ANOVA
experimental design and returns a pseudo-F value and a p-value. Since PERMANOVA
is nonparametric, significance is determined through permutations.

PERMDISP - This method analyzes the multivariate homogeneity of group
dispersions (variances). In essence, it determines whether the variances of
groups of samples are significantly different. The results of both parametric
and nonparametric significance tests are provided in the output. This method is
generally used as a companion to PERMANOVA. See vegan::betadisper for more
details.

db-RDA - This method is very similar to adonis and will only differ if certain
non-Euclidean semi- or non-metrics are used to generate the input distance
matrix, and negative eigenvalues are encountered. The only difference then will
be in the p-values, not the R^2 values. As part of the output, an ordination
plot is also generated that shows grouping/clustering of samples based on a
category in the metadata mapping file. This category is used to explain the
variability between samples. Thus, the ordination output of db-RDA is similar
to PCoA except that it is constrained, while PCoA is unconstrained (i.e. with
db-RDA, you must specify which category should be used to explain the
variability in your data). See vegan::capscale for more details.

For more information and examples pertaining to this script, please refer to
the accompanying tutorial, which can be found at
http://qiime.org/tutorials/category_comparison.html.
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("adonis example",
                                    "Runs the adonis statistical method on a distance matrix and mapping file "
                                    "using the Treatment category and 999 permutations, writing the output to the "
                                    "'adonis_out' directory.",
                                    "%prog --method adonis -i unweighted_unifrac_dm.txt -m Fasting_Map.txt -c "
                                    "Treatment -o adonis_out -n 999"))

script_info['script_usage'].append(("ANOSIM example",
                                    "Runs the ANOSIM statistical method on a distance matrix and mapping file "
                                    "using the Treatment category and 99 permutations, writing the output to the "
                                    "'anosim_out' directory.",
                                    "%prog --method anosim -i unweighted_unifrac_dm.txt -m Fasting_Map.txt -c "
                                    "Treatment -o anosim_out -n 99"))

script_info['output_description'] = """
At least one file will be created in the output directory specified by -o. For
most methods, a single output file containing the results of the test (e.g. the
effect size statistic and p-value) will be created. The format of the output
files will vary between methods as some are generated by native QIIME code,
while others are generated by R's vegan or ape packages. Please refer to the
script description for details on how to access additional information for
these methods, including what information is included in the output files.

db-RDA is the only exception in that two output files are created: a results
text file and a PDF of the ordination plot.
"""

script_info['required_options'] = [
    make_option('--method', help='the statistical method to use. Valid '
                'options: %s' % ', '.join(methods), type='choice', choices=methods),
    make_option('-i', '--input_dm', type='existing_filepath',
                help='the input distance matrix. WARNING: Only symmetric, hollow '
                'distance matrices may be used as input. Asymmetric distance '
                'matrices, such as those obtained by the UniFrac Gain metric (i.e. '
                'beta_diversity.py -m unifrac_g), should not be used as input'),
    make_option('-m', '--mapping_file', type='existing_filepath',
                help='the metadata mapping file'),
    make_option('-c', '--categories', type='string',
                help='a comma-delimited list of categories from the mapping file. '
                'Note: all methods except for BIO-ENV accept just a single category. If '
                'multiple categories are provided, only the first will be used'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    make_option('-n', '--num_permutations', help='the number of permutations '
                'to use when calculating statistical significance. Only applies to '
                'adonis, ANOSIM, MRPP, PERMANOVA, PERMDISP, and db-RDA. Must be '
                'greater than or equal to zero [default: %default]', default=999,
                type='int')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    out_dir = opts.output_dir
    categories = opts.categories.split(',')

    # Create the output dir if it doesn't already exist.
    try:
        if not exists(out_dir):
            create_dir(out_dir)
    except:
        option_parser.error("Could not create or access output directory '%s' "
                            "specified with the -o option." % out_dir)

    compare_categories(opts.input_dm, opts.mapping_file, opts.method,
                       categories, opts.num_permutations, out_dir)


if __name__ == "__main__":
    main()
