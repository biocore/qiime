#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Doug Wendel"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone", "Jesse Stombaugh", "Doug Wendel", \
                "Dan Knights", "Greg Caporaso", "Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Doug Wendel"
__email__ = "wendel@colorado.edu"
__status__ = "Release"
 

from os.path import isdir, isfile, join, basename, exists
from os import listdir, makedirs, getcwd
from glob import glob
from qiime.otu_category_significance import test_wrapper, test_wrapper_multiple,\
    sync_mapping_to_otu_table
from qiime.longitudinal_otu_category_significance import \
    longitudinal_otu_table_conversion_wrapper
from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.format import format_biom_table
from biom.parse import parse_biom_table, parse_biom_table_str
from qiime.parse import parse_mapping_file
from cogent.util.misc import create_dir

script_info={}
script_info['brief_description']="""OTU significance and co-occurence analysis"""
script_info['script_description']="""The script otu_category_significance.py  
tests whether any of the OTUs in an OTU table are significantly associated with 
a category in the category mapping file. This code uses ANOVA, the G test of 
independence, Pearson correlation, or a paired t-test to find OTUs that are 
differentially represented across experimental treatments or measured variables.

The script can also be used to measure co-occurrence. For instance it can also 
be used with presence/absence or abundance data for a phylogenetic group (such 
as that determined with quantitative PCR) to determine if any OTUs co-occur
with a taxon of interest, using the ANOVA, G test of Independence, or 
correlation.

One useful feature is to be able to run otu category significance across 
all taxonomic levels of an OTU table. For example, you can take your otu_table.biom 
and run summarize_taxa.py on it, which will create an output directory that 
contains your OTU table summarized at different taxonomic levels, with the 
file names containing L2, L3, L4 etc. to designate different taxonomic 
resolutions. This script allows you to "sweep" over those taxonomic levels 
contained within the directory by passing -i input_directory. This 
output will then contain an individual results file that 
is written for each summarized taxa table. Note for QIIME 1.6: the input  
directory that you pass to this script must contain biom tables. Thus, you 
might need to convert your summarized taxa tables that are in the classic OTU 
table form to biom tables using the convert_biom.py script. 

If you run multiple_rarefactions_even_depth.py on your OTU table, you may 
want to run otu_category_significance on each table and then collate the results 
achieved from each table. In order to do this, simply pass the directory as the 
input filepath, and then pass the -w collate option. See usage examples below. 

The statistical test to be run is designated with the -s option, and includes 
the following options:

The G test of independence (g_test): determines whether OTU presence/absence is
associated with a category (e.g. if an OTU is more or less likely to be present 
in samples from people with a disease vs healthy controls).

ANOVA (ANOVA): determines whether OTU relative abundance is different between 
categories (e.g. if any OTUs are increased or decreased in relative abundance in 
the gut microbiota of obese versus lean individuals). 

Pearson correlation (correlation): determines whether OTU abundance is 
correlated with a continuous variable in the category mapping file. (e.g. which
 OTUs are positively or negatively correlated with measured pH across soil samples)


The tests also include options for longitudinal data (i.e. datasets in which 
multiple samples are collected from a single individual or site.) The composition
of microbes may differ substantially across samples for reasons that do not 
relate to a study treatment. For instance, a given OTU may not be in an 
individual or study site for historical reasons, and so cannot change as a
result of a treatment. The longitudinal tests thus ignore samples from individuals
in which a particular OTU has never been observed across samples. The category 
mapping file must have an "individual" column indicating which sample is from 
which individual or site, and a "reference_sample" column, indicating which 
sample is the reference sample for an individual or site (e.g. time point zero 
in a timeseries experiment). The longitudinal options include:

Pearson correlation (longitudinal_correlation): determines whether OTU 
relative abundance is correlated with a continuous variable in the category 
mapping file while accounting for an experimental design where multiple samples
are collected from the same individual or site. Uses the change in relative
abundance for each sample from the reference sample (e.g. timepoint zero in 
a timeseries analysis) rather than the absolute relative abundances in the 
correlation (e.g. if the relative abundance before the treatment was 0.2, and 
after the treatment was 0.4, the new values for the OTU relative abundance will 
be 0.0 for the before sample, and 0.2 for the after, thus indicating that the 
OTU went up in response to the treatment.)

Paired t-test (paired_T): This option is when measurements were taken "before" 
and "after" a treatment. There must be exactly two measurements for each 
individual/site. The category mapping file must again have an individual column, 
indicating which sample is from which individual, and a reference_sample column, 
that has a 1 for the before time point and a 0 for the after.

With the exception of longitudinal correlation and paired_T, this script can be 
performed on a directory of OTU tables (for example, the output of 
multiple_rarefactions_even_depth.py), in addition to on a single OTU table. 
If the script is called on a directory, the resulting p-values are the average 
of the p-values observed when running a single test on each otu_table 
separately. It is generally a good practice to rarefy the OTU table 
(e.g. with single_rarefaction.py) prior to running these significance tests 
in order to avoid artifacts or biases from unequal sample sizes.
"""
script_info['script_usage']=[]

script_info['script_usage'].append(("G-test","""Perform a G test on otu_table.biom testing OTUs for differences in the abundance across the category "Treatment":""","""%prog -i otu_table.biom -m Fasting_Map.txt -s g_test -c Treatment -o single_g_test.txt"""))

script_info['script_usage'].append(("ANOVA","""Perform an ANOVA on otu_table.biom testing OTUs for differences in the abundance across the category "Treatment":""","""%prog -i otu_table.biom -m Fasting_Map.txt -s ANOVA -c Treatment -o single_anova.txt"""))

script_info['script_usage'].append(("ANOVA on multiple OTU tables and collate results","""Perform an ANOVA on all OTU tables in rarefied_otu_tables testing OTUs for differences in the abundance across the category "Treatment" and collate the results into one file:""","""%prog -i rarefied_otu_tables/ -m Fasting_Map.txt -s ANOVA -c Treatment -o multiple_anova.txt -w """))

script_info['script_usage'].append(("ANOVA on multiple OTU tables and write out separate results files","""Perform an ANOVA on all OTU tables in rarefied_otu_tables testing OTUs for differences in the abundance across the category "Treatment" and produce a results file for each OTU table:""","""%prog -i rarefied_otu_tables/ -m Fasting_Map.txt -s ANOVA -c Treatment -o rarefied_table_results/ """))

# do we have good input data for this? 
# script_info['script_usage'].append(("Example 2","""If the user would like to perform the same test using numerical qPCR data, where everything below a threshold value should be considered "absent" and everything above that value "present", the user will need to set the threshold by running the following command:""","""%prog -i otu_table.biom -m Mapping_file.txt -s g_test -c Treatment -t 0.16"""))

script_info['output_description']="""The G test results are output as tab delimited text, which can be examined in Excel. The output has the following columns:

* OTU: The name of the OTU.
* g_val: The raw test statistic.
* g_prob: The probability that this OTU is non-randomly distributed across the
 categories.
* Bonferroni_corrected: The probability after correction for multiple 
comparisons with the Bonferroni correction. In this correction, the p-value 
is multiplied by the number of comparisons performed (the number of OTUs 
    remaining after applying the filter).
* FDR_corrected: The probability after correction with the "false discovery rate"
 method. In this method, the raw p-values are ranked from low to high. Each 
 p-value is multiplied by the number of comparisons divided by the rank. This 
 correction is less conservative than the Bonferroni correction. The list of 
 significant OTUs is expected to have the percent of false positives predicted 
 by the p value.
* Contingency table columns: The next columns give the information in the 
contingency table and will vary in number and name based on the number of 
categories and their names. The two numbers in brackets represent the number of 
samples that were observed in those categories and the number that would be 
expected if the OTU members were randomly distributed across samples in the 
different categories. These columns can be used to evaluate the nature of a 
non-random association (e.g. if that OTU is always present in a particular 
    category or if it is never present).
* Consensus lineage: The consensus lineage for that OTU will be listed in the 
last column if it was present in the input OTU table.

The ANOVA results are output as tab delimited text that can be examined in 
Excel. The output has the following columns:

* OTU: The name of the OTU.
* prob: The raw probability from the ANOVA 
* Bonferroni_corrected: The probability after correction for multiple 
comparisons with the Bonferroni correction. In this correction, the p-value 
is multiplied by the number of comparisons performed (the number of OTUs 
    remaining after applying the filter). 
* FDR_corrected: The probability after correction with the "false discovery rate" 
method. In this method, the raw p-values are ranked from low to high. Each 
p-value is multiplied by the number of comparisons divided by the rank. This 
correction is less conservative than the Bonferroni correction. The list of 
significant OTUs is expected to have the percent of false positives predicted 
by the p value.
* Category Mean Columns: Contains one column for each category reporting the 
mean count of the OTU in that category.
* Consensus lineage: The consensus lineage for that OTU will be listed in the 
last column if it was present in the input OTU table.

The correlation and longitudinal_correlation test results are output as tab 
delimited text, which can be examined in Excel. The output has the following columns:

* OTU: The name of the OTU.  
* prob: The probability that the OTU relative abundance is correlated with the
 category values across samples. 
* otu_values_y: a list of the values (relative abundance) of the OTU across the
 samples that were plotted on the y axis for the correlation.
* cat_values_x: a list of the values of the selected category that were plotted
 on the x axis for the correlation.
* Bonferroni_corrected: The probability after correction for multiple
 comparisons with the Bonferroni correction. In this correction, the p-value
  is multiplied by the number of comparisons performed (the number of OTUs
   remaining after applying the filter). 
* FDR_corrected: The probability after correction with the "false discovery rate"
 method. In this method, the raw p-values are ranked from low to high. Each 
 p-value is multiplied by the number of comparisons divided by the rank. This 
 correction is less conservative than the Bonferroni correction. The list of 
 significant OTUs is expected to have the percent of false positives predicted
  by the p value.
* r: Pearson's r. This value ranges from -1 to +1, with -1 indicating a perfect
 negative correlation, +1 indicating a perfect positive correlation, and 0 
 indicating no relationship.
* Consensus lineage: The consensus lineage for that OTU will be listed in the 
last column if it was present in the input OTU table.

The paired_T results are output as tab delimited text that can be examined in
 Excel. The output has the following columns:

* OTU: The name of the OTU.
* prob: The raw probability from the paired T test
* T stat: The raw T value
* average_diff: The average difference between the before and after samples in 
the individuals in which the OTU was observed.
* num_pairs: The number of sample pairs (individuals) in which the OTU was observed.
* Bonferroni_corrected: The probability after correction for multiple comparisons
 with the Bonferroni correction. In this correction, the p-value is multiplied 
 by the number of comparisons performed (the number of OTUs remaining after 
    applying the filter). 
* FDR_corrected: The probability after correction with the "false discovery 
rate" method. In this method, the raw p-values are ranked from low to high. 
Each p-value is multiplied by the number of comparisons divided by the rank. 
This correction is less conservative than the Bonferroni correction. The list 
of significant OTUs is expected to have the percent of false positives 
predicted by the p value.
* Consensus lineage: The consensus lineage for that OTU will be listed in the
 last column if it was present in the input OTU table.

"""

script_info['required_options']=[
    make_option('-i','--otu_table_fp',
        help='path to the otu table in biom format, or to a directory '  
             'containing OTU tables',type='existing_path'),
    make_option('-m','--category_mapping_fp',type='existing_filepath',
        help='path to category mapping file'),
    make_option('-o', '--output_fp', dest='output_fp', help='path to the output '
                'file or directory', type='new_filepath')
]

script_info['optional_options']=[
    make_option('-c', '--category', dest='category', type='string', default=None,
        help='name of the category over which to run the analysis'),
    make_option('-s','--test', dest='test', default='ANOVA',
        help='the type of statistical test to run. options are: '
        'g_test: determines whether OTU presence/absence is associated '
        'with a category using the G test of Independence.      '
        'ANOVA: determines whether OTU abundance is associated with a '
        'category.      ' 
        'correlation: determines whether OTU abundance is correlated ' 
        'with a continuous variable in the category mapping file.     ' 
        'longitudinal_correlation: determine whether OTU relative ' 
        'abundance is correlated with a continuous variable in the ' 
        'category mapping file in longitudinal study designs such as ' 
        'with timeseries data.     paired_T: determine whether OTU ' 
        'relative abundance goes up or down in response to a treatment. [default: %default]',
        type="choice",choices=["g_test", "ANOVA", "correlation", 
        "longitudinal_correlation", "paired_T"]),
    make_option('-f','--filter', dest='filter', type='float', 
        default= 0.25, 
        help='minimum fraction of samples that must contain the OTU for the ' 
        'OTU to be included in the analysis. For longitudinal options, is ' 
        'the fraction of individuals/sites that were not ' 
        'ignored because of the OTU not being observed in any of the ' 
        'samples from that individual/site. [default: %default]'),
    make_option('-t','--threshold', dest='threshold', default=None, type='float', 
        help='threshold under which to consider something absent: ' 
        'Only used if you have numerical data that should be converted to ' 
        'present or absent based on a threshold. Should be None for ' 
        'categorical data or with the correlation test. default value is None'),
    make_option('-l', '--otu_include_fp', type='existing_filepath', default=None,
        help='path to a file with a list of OTUs to evaluate. By default ' 
        'evaluates all OTUs that pass the minimum sample filter. If a ' 
        'filepath is given here in which each OTU name one wishes to ' 
        'evaluate is on a separate line, will apply this additional filter'),
    make_option('-z', '--reference_sample_column', dest='reference_sample_column',
        default='reference_sample', type='string',
        help='This column specifies the sample to which all other samples ' 
        'within an individual are compared. For instance, for timeseries ' 
        'data, it would usually be the initial timepoint before a treatment ' 
        'began. The reference samples should be marked with a 1, and ' 
        'other samples with a 0.'),
    make_option('-n','--individual_column', dest='individual_column', 
        default='individual', type='string',
        help='name of the column in the category mapping file that ' 
        'designates which sample is from which individual.'),
    make_option('-b', '--converted_otu_table_output_fp', type='new_filepath',
        default=None, help='the test options longitudinal_correlation and ' 
        'paired_T convert the original OTU table into one in which samples ' 
        'that are ignored because they are never observed in an individual ' 
        'are replaced with the ignore number 999999999 and the OTU counts are ' 
        'the change in relative abundance compared to the designated reference ' 
        'sample. If a filepath is given with the -b option ' 
        'this converted OTU table will be saved to this path.'),
    make_option('--relative_abundance', default=False, help='Some of the '
        'statistical tests, such as Pearson correlation and ANOVA, convert '
        'the OTU counts to relative abundances prior to performing the '
        'calculations. This parameter can be set if a user wishes to disable '
        'this step. (e.g. if an OTU table has already been converted '
        'to relative abundances.)'),
    make_option('-w', '--collate_results', dest='collate_results',
        action='store_true', default=False,
        help='When passing in a directory of OTU tables, '
        'this parameter gives you the option of collating those resulting '
        'values. For example, if your input directory contained multiple '
        'rarefied OTU tables at the same depth, pass the -w option '
        'in order to find the average p-value for your statistical test '
        'over all rarefied tables and collate the results into one file. '
        'If your input directory contained OTU tables '
        'that contained different taxonomic levels, filtering levels, etc '
        'then do not pass the -w option so that an individual results file '
        'is created for every input OTU table. [default: %default]')]
    
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    otu_table_fp = opts.otu_table_fp
    otu_include_fp = opts.otu_include_fp
    output_fp = opts.output_fp
    verbose = opts.verbose
    category_mapping_fp = opts.category_mapping_fp
    individual_column = opts.individual_column
    reference_sample_column = opts.reference_sample_column
    conv_output_fp = opts.converted_otu_table_output_fp
    relative_abundance = opts.relative_abundance
    filter = opts.filter
    test = opts.test
    category = opts.category
    threshold = opts.threshold
    collate_results = opts.collate_results

    # check mapping file
    category_mapping = open(category_mapping_fp,'U')
    category_mapping = parse_mapping_file(category_mapping)
    if not category:
        if test != 'paired_T':
            option_parser.error('a category in the category mapping file must be'
                ' specified with the -c option for this test')

    # set up threshold value for filtering, if any
    if threshold and threshold != 'None':
        threshold = float(threshold)

    # if specifying a list of OTUs to look at specifically
    if otu_include_fp and otu_include_fp != 'None':
        otu_include = open(otu_include_fp)
    else:
        otu_include = None

    # if only passing in a single OTU table
    if isdir(otu_table_fp) is False:
        # raise error if collate option is being passed to single table
        if collate_results is True:
            option_parser.error('Cannot collate the results of only one table.'
                ' Please rerun the command without passing the -w option')
        else:
            #open and parse the biom table fp
            otu_table = parse_biom_table(open(otu_table_fp, 'U'))

            # run the statistical test
            output = test_wrapper(test, otu_table, category_mapping, 
                            category, threshold, filter, otu_include, 
                            otu_table_relative_abundance=relative_abundance)
            
            # write output
            output_file = open(output_fp, 'w')
            output_file.write('\n'.join(output))
            output_file.close()

    # if the user has passed in a directory
    if isdir(otu_table_fp) is True:

        # negate_collate to return an results file on a per table basis
        if collate_results is False:
      
            # build list of otu tables
            otu_table_paths = glob('%s/*biom' % otu_table_fp)

            # if output dir doesn't exist, then make it
            if exists(output_fp):
                pass
            else:
                makedirs(output_fp)

            for otu_table_fp in otu_table_paths:
                #open and parse the biom table fp
                otu_table = parse_biom_table(open(otu_table_fp, 'U'))

                #synchronize the mapping file with the otu table
                category_mapping, removed_samples = \
                    sync_mapping_to_otu_table(otu_table, category_mapping)
                if removed_samples:
                    print "Warning, the following samples were in the category mapping file " +\
                                    "but not the OTU table and will be ignored: "
                    for i in removed_samples:
                            print i + '\n'

                # create naming convention for output file
                # will look like: otu_table_ANOVA_Treatment.txt
                output_basename = basename(otu_table_fp)
                output_basename = output_basename.replace(".biom","")
                output_fp_sweep = "%s_%s_%s.txt" % \
                    (output_basename,test,category)

                # if the convert_otu_table_fp is passed, save the converted table
                if test == 'longitudinal_correlation' or test == 'paired_T':
                    converted_otu_table = longitudinal_otu_table_conversion_wrapper(table, 
                        category_mapping, individual_column, reference_sample_column)
                    if conv_output_fp:
                        of = open(conv_output_fp, 'w')
                        of.write(format_biom_table(converted_otu_table))
                        of.close()
                    if test == 'longitudinal_correlation':
                        #set the otu_include list to all of the OTUs, this effectively
                        #deactivates the filter for correlation, because the filtered OTU_list is
                        #rewritten with the otu_include list in the test_wrapper
                        if not otu_include:
                            otu_include = set(otu_table.ObservationIds)
                        output = test_wrapper('correlation', converted_otu_table, 
                            category_mapping, category, threshold, filter, otu_include, 
                            999999999.0, True)
                    elif test == 'paired_T':
                        output = test_wrapper('paired_T', converted_otu_table, 
                            category_mapping, category, threshold, 
                            filter, otu_include, 999999999.0, True, 
                            individual_column, reference_sample_column)

                # run test single input table from the directory  
                else:
                    output = test_wrapper(test, otu_table, category_mapping, 
                        category, threshold, filter, otu_include, 
                        otu_table_relative_abundance=relative_abundance)
                    
                # write output file with new naming convention
                output_file = open(join(output_fp,output_fp_sweep), 'w')
                output_file.write('\n'.join(output))
                output_file.close()


        # Use when the input dir contains rarefied OTU tables, and you want
        # to collate the p-values & results into one results file       
        if collate_results is True:
            if test != 'longitudinal_correlation' and test != 'paired_T':
                # get biom tables
                otu_table_paths = glob('%s/*biom' % otu_table_fp)

                #get aggregated tables
                parsed_otu_tables= []
                for otu_table_fp in otu_table_paths:
                    otu_table = open(otu_table_fp, 'U')
                    otu_table = parse_biom_table(otu_table)
                    parsed_otu_tables.append(otu_table)

                #synchronize the mapping file with the otu table
                #checks with just the first OTU table and assumes that all otu tables
                #have the same collection of samples
                category_mapping, removed_samples = \
                    sync_mapping_to_otu_table(parsed_otu_tables[0],category_mapping)
                if removed_samples:
                    print "Warning, the following samples were in the category mapping file " +\
                                "but not the OTU table and will be ignored: "
                    for i in removed_samples:
                            print i + '\n'

                # get output from statistical test            
                output = test_wrapper_multiple(test, parsed_otu_tables, 
                    category_mapping, category, threshold, filter, otu_include,
                    otu_table_relative_abundance=relative_abundance)

                #write out aggregated results
                output_file = open(output_fp, 'w')
                output_file.write('\n'.join(output))
                output_file.close()

            else:
                option_parser.error("You cannot collate the results obtained from "
                    "using the longitudinal_correlation and paired_T options.")

if __name__ == "__main__":
    main()