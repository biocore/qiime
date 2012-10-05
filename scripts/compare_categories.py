#!/usr/bin/env python
from __future__ import division

__author__ = "Logan Knecht"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Logan Knecht", "Michael Dwan", "Damien Coy", "Jai Ram Rideout",
               "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from os import path, makedirs, listdir

from cogent.util.misc import create_dir
from numpy import zeros
from numpy.random import permutation

from qiime.format import format_p_value_for_num_iters
from qiime.parse import (parse_distmat, fields_to_dict, parse_mapping_file,
                         parse_mapping_file_to_dict)
from qiime.stats import Anosim, BioEnv, Permanova
from qiime.util import (parse_command_line_parameters, make_option,
                        get_options_lookup, DistanceMatrix, MetadataMap,
                        RExecutor)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Analyzes statistical significance of sample groupings using distance matrices"""
script_info['script_description'] = """
This script allows for the analysis of the strength and statistical
significance of sample groupings using a distance matrix as the primary input.
Several statistical methods are available: adonis, ANOSIM, BEST, Moran's I,
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

BEST - This method looks at the numerical environmental variables relating
samples in a distance matrix. For instance, given a UniFrac distance matrix and
pH and latitude (or any other number of variables) in soil samples, BEST will
rank them in order of which best explain patterns in the communities. This
method will only accept categories that are numerical (continuous or discrete).
This is currently the only method in this script that accepts more than one
category (via -c). See vegan::bioenv for more details.

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
script_info['script_usage'].append(("adonis",
"Runs the adonis statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category and 999 permutations. Then it outputs the "
"results to the 'adonis' directory. The full file path will be: "
"./adonis/adonis_results.txt",
"%prog --method adonis -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o adonis -n 999"))

script_info['script_usage'].append(("ANOSIM",
"Runs the ANOSIM statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category and 999 perutations. Then it outputs the "
"results to the 'anosim' directory. The full file path will be: "
"./anosim/anosim_results.txt",
"%prog --method anosim -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim -n 999"))

script_info['script_usage'].append(("BEST",
"Runs the BEST statistical method on a distance matrix and mapping file "
"using the LATITUDE and LONGITUDE categories. Then it outputs the results to "
"the 'best' directory. The full file path will be: ./best/best_results.txt",
"%prog --method best -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c LATITUDE,LONGITUDE -o best"))

script_info['script_usage'].append(("Moran's I",
"Runs the Moran's I statistical method on a distance matrix and mapping "
"file using the PH category. Then it outputs the results to the 'morans_i' "
"directory. The full file path will be: ./morans_i/Morans_I_results.txt",
"%prog --method morans_i -i  datasets/88_soils/unweighted_unifrac_dm.txt -m "
"datasets/88_soils/map.txt -c PH -o morans_i"))

script_info['script_usage'].append(("MRPP",
"Runs the MRPP statistical method on a distance matrix and mapping file using "
"the HOST_SUBJECT_ID category. Then it outputs the results to the 'mrpp' "
"directory. The full file path will be: ./mrpp/mrpp_results.txt",
"%prog --method mrpp -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o mrpp -n 999"))

script_info['script_usage'].append(("PERMANOVA", "Runs the PERMANOVA "
"statistical method on a distance matrix and mapping file using the "
"HOST_SUBJECT_ID category. Then it outputs the results to the 'permanova' "
"directory. The full file path will be: ./permanova/permanova_results.txt",
"%prog --method permanova -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova -n 999"))

script_info['script_usage'].append(("PERMDISP", "Runs the PERMDISP "
"statistical method on a distance matrix and mapping file using the "
"HOST_SUBJECT_ID category. Then it outputs the results to the 'permdisp' "
"directory. The full file path will be: ./permdisp/permdisp_results.txt",
"%prog --method permdisp -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permdisp"))

script_info['script_usage'].append(("db-RDA", "Runs the db-RDA statistical "
"method on a distance matrix and mapping file using the HOST_SUBJECT_ID "
"category. Then it outputs the results to the 'dbrda' directory. The full "
"file path will be: ./dbrda/dbrda_results.txt and ./dbrda/dbrda_plot.txt",
"%prog --method dbrda -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o dbrda"))

script_info['output_description']= """
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
        'options: [adonis, anosim, best, morans_i, mrpp, permanova, '
        'permdisp, dbrda]', type='choice', choices=['adonis', 'anosim', 'best',
        'morans_i', 'mrpp', 'permanova', 'permdisp', 'dbrda']),
    make_option('-i', '--input_dm', type='existing_filepath',
        help='the input distance matrix. WARNING: Only symmetric, hollow '
        'distance matrices may be used as input. Asymmetric distance '
        'matrices, such as those obtained by the UniFrac Gain metric (i.e. '
        'beta_diversity.py -m unifrac_g), should not be used as input'),
    make_option('-m', '--mapping_file', type='existing_filepath',
        help='the metadata mapping file'),
    make_option('-c', '--categories', type='string',
        help='a comma-delimited list of categories from the mapping file. '
        'Note: all methods except for BEST take just a single category. If '
        'multiple categories are provided, only the first will be used'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    make_option('-n', '--num_permutations', help='the number of permutations '
        'to use when calculating statistical significance. Only applies to '
        'adonis, ANOSIM, MRPP, PERMANOVA, PERMDISP, and db-RDA '
        '[default: %default]', default=999, type='int')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    try:
        if not path.exists(opts.output_dir):
            create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    # Parse the mapping file and distance matrix.
    md_map = MetadataMap.parseMetadataMap(open(opts.mapping_file,'U'))
    dm = DistanceMatrix.parseDistanceMatrix(open(opts.input_dm,'U'))

    # Separate all categories into a list, then grab the first category.
    categories = opts.categories.split(',')

    # Cursory check to make sure all categories passed in are in mapping file.
    maps = parse_mapping_file(open(opts.mapping_file,'U').readlines())
    for category in categories:
        if not category in maps[1][1:]:
            option_parser.error("Category '%s' not found in mapping file "
                                "columns:" % category)

    # Make sure the input distance matrix is symmetric and hollow. Must check
    # here before allowing R to use it, as R will silently ignore the diagonal
    # and upper triangle of the distance matrix.
    if not dm.is_symmetric_and_hollow():
        option_parser.error("The distance matrix must be symmetric and "
                            "hollow.")

    # Figure out which method we need to run.
    if opts.method == 'adonis':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + categories[0] + " -o " + opts.output_dir + " -n " + \
            str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "adonis.r", output_dir=opts.output_dir)
    elif opts.method == 'anosim':
        anosim = Anosim(md_map, dm, categories[0])
        anosim_results = anosim(opts.num_permutations)

        output_file = open(opts.output_dir + "/" + opts.method + \
            "_results.txt", "w+")
        output_file.write("Method Name\tR-value\tP-value")
        output_file.write("\n")
        output_file.write(anosim_results["method_name"]+"\t"+\
            str(anosim_results["r_value"])+"\t"+\
            str(anosim_results["p_value"])+"\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'best':
        bioenv = BioEnv(dm, md_map, categories)
        bioenv_results = bioenv()

        output_file = open(opts.output_dir+"/best_results.txt", 'w+')
        output_file.write("Method Name:\tNum_Vars:\t")
        output_file.write("\n")
        output_file.write(bioenv_results["method_name"]+"\t"+\
            str(bioenv_results["num_vars"]) + "\t")
        output_file.write("\n")
        output_file.write("Variables:\t")
        output_file.write("\n")
        for variable in bioenv_results["vars"]:
            output_file.write(str(variable) + "\t")
        output_file.write("\n")
        output_file.write("RHO_Values:\t")
        output_file.write("\n")
        for rho_val in bioenv_results["bioenv_rho_vals"]:
            output_file.write(str(rho_val) + "\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'morans_i':
        command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + categories[0] + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "morans_i.r", output_dir=opts.output_dir)
    elif opts.method == 'mrpp':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + categories[0] + " -o " + opts.output_dir + \
            " -n " + str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "mrpp.r", output_dir=opts.output_dir)
    elif opts.method == 'permanova':
        permanova_plain = Permanova(md_map, dm, categories[0])
        permanova_results = permanova_plain(opts.num_permutations)

        output_file = open(opts.output_dir+"/permanova_results.txt", 'w+')
        output_file.write("Method Name\tF-value\tP-value")
        output_file.write("\n")
        output_file.write(permanova_results["method_name"]+"\t"+\
            str(permanova_results["f_value"]) + "\t" + \
            format_p_value_for_num_iters(permanova_results["p_value"], \
            opts.num_permutations)+"\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'permdisp':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + categories[0] + " -o " + opts.output_dir + " -n " + \
            str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "permdisp.r", output_dir=opts.output_dir)
    elif opts.method == 'dbrda':
        command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + categories[0] + " -o " + opts.output_dir + " -n " + \
            str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "dbrda.r", output_dir=opts.output_dir)


if __name__ == "__main__":
    main()
