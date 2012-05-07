#!/usr/bin/env python
from __future__ import division

__author__ = "Logan Knecht"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Logan Knecht", "Michael Dwan", "Damien Coy", "Jai Ram Rideout",
               "Levi McCracken"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
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
script_info['brief_description'] = """
Analyzes distance matrices for statistical significance of sample grouping
"""
script_info['script_description'] = """
This script allows for the analysis of distance matrices using several \
statistical methods. These methods are Adonis, Anosim, BEST, DFA, Moran's I, \
MRPP, PERMANOVA, PERMDISP, RDA.

Adonis - This method takes a distance matrix and mapping file. It then \
identifies important points in the data and performs F-tests on the initial \
data, and random permutations (via shuffling) the category data. Then, \
it finally returns the information that was identified in the samples. It's \
stated that it partitions (or seperates the data) for this analysis in order\
 to find underlying relationships.

Anosim - This method takes in a distance matrix and a mapping file. \
This method tests whether two or more categories are significantly \
different. You can specify a category in the metadata mapping file to separate \
samples into groups and then test whether there are significant differences \
between those groups. For example, you might test whether Control samples are \
significantly different from Fast samples. Since ANOSIM is non-parametric, \
significance is determined through permutations.

BEST - This method looks at the numerical environmental variables relating \
samples in a distance matrix. For instance, the unifrac distance matrix and \
pH and latitude (or any other number of variables) in soil samples, and ranks \
them in order of which best explain patterns in the communities. This method \
will only accept categories that are numerical.

Moran's I - This method takes in a distance matrix and mapping file. Then it \
uses the geographical data supplied to identify what type of spatial \
configuration occurs in the samples. Are they dispersed, clustered, or of no \
distinctly noticeable configuration when compared to each other? This method \
will only accept a category that is numerical.

MRPP - This method takes in a distance matrix and a mapping file. It then \
tests whether two or more categories are significantly different. You can \
specify a category in the metadata mapping file to separate samples into \
groups and then test whether there are significant differences between those \
groups. For example, you might test whether Control samples are significantly \
different from Fast samples. Since MRPP is non-parametric, significance is \
determined through permutations.

PERMANOVA - This method takes distance matrix and a mapping file. This method \
is for testing the simultaneous response of one or more variables to one or \
more factors in an ANOVA experimental design on the basis of any distance \
metric. It returns a R value and a P value.

PERMDISP - This method takes a distance matrix and a mapping file. \
This method is a procedure for the analysis of multivariate homogeneity of \
group dispersions (variances). Permutations can be utilized to measure the \
dissimilatities between groups.

RDA - This method takes a distance matrix and a mapping file. This method \
is an ordination method that shows grouping/clustering of samples based on \
a category in the metadata mapping file and a distance matrix. This category \
is used to explain the variability between samples. Thus, RDA is similar to \
PCoA except that it is constrained, while PCoA is unconstrained (you must \
specify which category should be used to explain the variability in your data).

For more information and examples pertaining to this script, please refer to \
the accompanying tutorial, which can be found at \
http://qiime.org/tutorials/category_comparison.html.
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("Adonis",
"Performs the Adonis statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category and 999 permutations. Then it outputs the "
"results to the 'adonis' directory. The full file path will be: "
"./adonis/adonis_results.txt",
"%prog --method adonis -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o adonis -n 999"))

script_info['script_usage'].append(("Anosim",
"Performs the Anosim statistical method on a distance matrix and mapping file "
"using the HOST_SUBJECT_ID category and 999 perutations. Then it outputs the "
"results to the 'anosim' directory. The full file path will be: "
"./anosim/anosim_results.txt",
"%prog --method anosim -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o anosim -n 999"))

script_info['script_usage'].append(("BEST",
"Performs the BEST statistical method on a distance matrix and mapping file "
"using the LATITUDE and LONGITUDE categories. Then it outputs the results to "
"the 'best' directory. The full file path will be: ./best/best_results.txt",
"%prog --method best -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c LATITUDE,LONGITUDE -o best"))

script_info['script_usage'].append(("Moran's I",
"Performs the Moran's I statistical method on a distance matrix and mapping "
"file using the PH category. Then it outputs the results to the 'morans_i' "
"directory. The full file path will be: ./morans_i/Morans_I_results.txt",
"%prog --method morans_i -i  datasets/88_soils/unweighted_unifrac_dm.txt -m "
"datasets/88_soils/map.txt -c PH -o morans_i"))

script_info['script_usage'].append(("MRPP", "Performs the MRPP statistical "
"method on a distance matrix and mapping file using the HOST_SUBJECT_ID "
"category. Then it outputs the results to the 'mrpp' directory. The full file "
"path will be: ./mrpp/mrpp_results.txt",
"%prog --method mrpp -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o mrpp -n 999"))

script_info['script_usage'].append(("PERMANOVA", "Performs the PERMANOVA "
"statistical method on a distance matrix and mapping file using the "
"HOST_SUBJECT_ID category. Then it outputs the results to the 'permanova' "
"directory. The full file path will be: ./permanova/permanova_results.txt",
"%prog --method permanova -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permanova -n 999"))

script_info['script_usage'].append(("PERMDISP", "Performs the PERMDISP "
"statistical method on a distance matrix and mapping file using the "
"HOST_SUBJECT_ID category. Then it outputs the results to the 'permdisp' "
"directory. The full file path will be: ./permdisp/betadisper_results.txt",
"%prog --method permdisp -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o permdisp"))

script_info['script_usage'].append(("RDA", "Performs the RDA statistical "
"method on a distance matrix and mapping file using the HOST_SUBJECT_ID "
"category. Then it outputs the results to the 'rda' directory. The full "
"file path will be: ./RDA/rda_results.txt and ./RDA/rda_plot.txt",
"%prog --method rda -i datasets/keyboard/unweighted_unifrac_dm.txt -m "
"datasets/keyboard/map.txt -c HOST_SUBJECT_ID -o rda"))

script_info['output_description']= """
Adonis:
One file is created and outputs the results into it. The results will be: \
Analysis of variance(AOV) table, degrees of freedom, sequential sums of \
squares, mean squares, F statistics, partial R-squared and p-values, based \
on the N permutations.

Anosim:
One file is output to the designated location under the name of \
anosim_results.txt. The information in the file will be an R-value and \
p-value.

Best:
This outputs one file 'best_results.txt' It will have teh method name, \
The number of variables. The list of varibles supplied. And lastly, the \
RHO values, which are ranked pearson correlations for the best combination \
of variables that describe the community.

Moran's I:
The output file is placed in a directory specified by -o. The file will be \
a text file with 4 values: observed, expected, sd, and p.value. The observed \
value is Morans I index of x. This is computed based on the values passed in \
to be compared with the weights. The expected value is the value of I under \
the null hypothesis. The sd is the standard deviation of I under the null \
hypothesis. P Value is the p-value of the test of the null hypothesis against \
the alternative hypothesis specified in alternative. Each of these values, \
except for the p-value, should be between -1 and 1.

MRPP:
The command in the previous section creates a single output file in the \
directory specified by the -o arguement, if not it will be sent to the \
directory location it was called from. The file will be named \
mrpp_results.txt. The file will contain a dissimilarity index, the class mean \
and counts. It will also conatin information about the chance corrected \
within-group agreement A, as well as the result Based on observed delta, and \
expected delta. There will also be the Significance of delta and the amount \
of permutations performed.

PERMANOVA:
Permanova returns one output file containing the the file passed in, the \
F-value and the p-value.

PERMDISP:
This method returns one file that outputs an analysis of varaiance table. \
Responses with the distances will be shown. There will be the strata \
relationship, then the sample information as well. Lastly you will be \
able to see the f-value and p-value.

RDA:
RDA outputs a two files. One is calles rda_results.txt, the other file \
is rda_plot.pdf. rda.txt contains the Inertia Proportion Rank, the \
Eigenvalues for constrained axes, and the Eigenvalues for unconstrained axes.
"""

script_info['required_options'] = [
    # All methods use these
    make_option('--method', help='The category analysis method. Valid '
        'options: [adonis, anosim, best, morans_i, mrpp, permanova, '
        'permdisp, rda]', type='choice', choices=['adonis', 'anosim', 'best',
        'morans_i', 'mrpp', 'permanova', 'permdisp', 'rda']),
    make_option('-i', '--input_dm', help='the input distance matrix'),
    make_option('-m', '--mapping_file', help='the metadata mapping file'),
    make_option('-c', '--categories', help='A comma-delimited list of '
        'categories from the mapping file (NOTE: many methods take just a '
        'single category, if multiple are passed only the first will be '
        'selected.)'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    # Only some methods use permutations.
    make_option('-n', '--num_permutations', help='the number of permutations '
        'to perform. Only applies to adonis, anosim, mrpp, and permanova '
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
    md_map = MetadataMap.parseMetadataMap(open(opts.mapping_file))
    dm = DistanceMatrix.parseDistanceMatrix(open(opts.input_dm))

    # Separate all categories into a list, then grab the first category.
    categories = opts.categories.split(',')
    first_category = categories[0]

    # Cursory check to make sure all categories passed in are in mapping file.
    maps = parse_mapping_file(open(opts.mapping_file,'U').readlines())
    for category in categories:
        if not category in maps[1][1:]:
            option_parser.error("Category '%s' not found in mapping file "
                                "columns:" % category)
    # Figure out which method we need to run.
    if opts.method == 'adonis':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir + " -n " + \
            str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "adonis.r", output_dir=opts.output_dir)
    elif opts.method == 'anosim':
        anosim = Anosim(md_map, dm, first_category)
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
            " -c " + first_category + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "morans_i.r", output_dir=opts.output_dir)
    elif opts.method == 'mrpp':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir + \
            " -n " + str(opts.num_permutations)]
        rex = RExecutor()
        rex(command_args, "mrpp.r", output_dir=opts.output_dir)
    elif opts.method == 'permanova':
        permanova_plain = Permanova(md_map, dm, first_category)
        permanova_results = permanova_plain(opts.num_permutations)

        output_file = open(opts.output_dir+"/permanova_results.txt", 'w+')
        output_file.write("Method Name\tR-value\tP-value")
        output_file.write("\n")
        output_file.write(permanova_results["method_name"]+"\t"+\
            str(permanova_results["r_value"]) + "\t" + \
            format_p_value_for_num_iters(permanova_results["p_value"], \
            opts.num_permutations)+"\t")
        output_file.write("\n")
        output_file.close()
    elif opts.method == 'permdisp':
        command_args = ["-d " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "betadisper.r", output_dir=opts.output_dir)
    elif opts.method == 'rda':
        command_args = ["-i " + opts.input_dm + " -m " + opts.mapping_file + \
            " -c " + first_category + " -o " + opts.output_dir]
        rex = RExecutor()
        rex(command_args, "rda.r", output_dir=opts.output_dir)

if __name__ == "__main__":
    main()
