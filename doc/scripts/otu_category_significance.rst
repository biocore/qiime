.. _otu_category_significance:

.. index:: otu_category_significance.py

*otu_category_significance.py* -- OTU significance and co-occurence analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The script `otu_category_significance.py <./otu_category_significance.html>`_  
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
and run `summarize_taxa.py <./summarize_taxa.html>`_ on it, which will create an output directory that 
contains your OTU table summarized at different taxonomic levels, with the 
file names containing L2, L3, L4 etc. to designate different taxonomic 
resolutions. This script allows you to "sweep" over those taxonomic levels 
contained within the directory by passing -i input_directory. This 
output will then contain an individual results file that 
is written for each summarized taxa table. Note for QIIME 1.6: the input  
directory that you pass to this script must contain biom tables. Thus, you 
might need to convert your summarized taxa tables that are in the classic OTU 
table form to biom tables using the `convert_biom.py <http://biom-format.org/documentation/biom_conversion.html>`_ script. 

If you run `multiple_rarefactions_even_depth.py <./multiple_rarefactions_even_depth.html>`_ on your OTU table, you may 
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
`multiple_rarefactions_even_depth.py <./multiple_rarefactions_even_depth.html>`_), in addition to on a single OTU table. 
If the script is called on a directory, the resulting p-values are the average 
of the p-values observed when running a single test on each otu_table 
separately. It is generally a good practice to rarefy the OTU table 
(e.g. with `single_rarefaction.py <./single_rarefaction.html>`_) prior to running these significance tests 
in order to avoid artifacts or biases from unequal sample sizes.



**Usage:** :file:`otu_category_significance.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the otu table in biom format, or to a directory containing OTU tables
	-m, `-`-category_mapping_fp
		Path to category mapping file
	-o, `-`-output_fp
		Path to the output file or directory
	
	**[OPTIONAL]**
		
	-c, `-`-category
		Name of the category over which to run the analysis
	-s, `-`-test
		The type of statistical test to run. options are: g_test: determines whether OTU presence/absence is associated with a category using the G test of Independence.      ANOVA: determines whether OTU abundance is associated with a category.      correlation: determines whether OTU abundance is correlated with a continuous variable in the category mapping file.     longitudinal_correlation: determine whether OTU relative abundance is correlated with a continuous variable in the category mapping file in longitudinal study designs such as with timeseries data.     paired_T: determine whether OTU relative abundance goes up or down in response to a treatment. [default: ANOVA]
	-f, `-`-filter
		Minimum fraction of samples that must contain the OTU for the OTU to be included in the analysis. For longitudinal options, is the fraction of individuals/sites that were not ignored because of the OTU not being observed in any of the samples from that individual/site. [default: 0.25]
	-t, `-`-threshold
		Threshold under which to consider something absent: Only used if you have numerical data that should be converted to present or absent based on a threshold. Should be None for categorical data or with the correlation test. default value is None
	-l, `-`-otu_include_fp
		Path to a file with a list of OTUs to evaluate. By default evaluates all OTUs that pass the minimum sample filter. If a filepath is given here in which each OTU name one wishes to evaluate is on a separate line, will apply this additional filter
	-z, `-`-reference_sample_column
		This column specifies the sample to which all other samples within an individual are compared. For instance, for timeseries data, it would usually be the initial timepoint before a treatment began. The reference samples should be marked with a 1, and other samples with a 0.
	-n, `-`-individual_column
		Name of the column in the category mapping file that designates which sample is from which individual.
	-b, `-`-converted_otu_table_output_fp
		The test options longitudinal_correlation and paired_T convert the original OTU table into one in which samples that are ignored because they are never observed in an individual are replaced with the ignore number 999999999 and the OTU counts are the change in relative abundance compared to the designated reference sample. If a filepath is given with the -b option this converted OTU table will be saved to this path.
	`-`-relative_abundance
		Some of the statistical tests, such as Pearson correlation and ANOVA, convert the OTU counts to relative abundances prior to performing the calculations. This parameter can be set if a user wishes to disable this step. (e.g. if an OTU table has already been converted to relative abundances.)
	-w, `-`-collate_results
		When passing in a directory of OTU tables, this parameter gives you the option of collating those resulting values. For example, if your input directory contained multiple rarefied OTU tables at the same depth, pass the -w option in order to find the average p-value for your statistical test over all rarefied tables and collate the results into one file. If your input directory contained OTU tables that contained different taxonomic levels, filtering levels, etc then do not pass the -w option so that an individual results file is created for every input OTU table. [default: False]


**Output:**

The G test results are output as tab delimited text, which can be examined in Excel. The output has the following columns:

* OTU: The name of the OTU.
* g_val: The raw test statistic.
* g_prob: The probability that this OTU is non-randomly distributed across the categories.
* Bonferroni_corrected: The probability after correction for multiple comparisons with the Bonferroni correction. In this correction, the p-value is multiplied by the number of comparisons performed (the number of OTUs remaining after applying the filter).
* FDR_corrected: The probability after correction with the "false discovery rate" method. In this method, the raw p-values are ranked from low to high. Each p-value is multiplied by the number of comparisons divided by the rank. This correction is less conservative than the Bonferroni correction. The list of significant OTUs is expected to have the percent of false positives predicted by the p value.
* Contingency table columns: The next columns give the information in the contingency table and will vary in number and name based on the number of categories and their names. The two numbers in brackets represent the number of samples that were observed in those categories and the number that would be expected if the OTU members were randomly distributed across samples in the different categories. These columns can be used to evaluate the nature of a non-random association (e.g. if that OTU is always present in a particular category or if it is never present).
* Consensus lineage: The consensus lineage for that OTU will be listed in the last column if it was present in the input OTU table.

The ANOVA results are output as tab delimited text that can be examined in Excel. The output has the following columns:

* OTU: The name of the OTU.
* prob: The raw probability from the ANOVA 
* Bonferroni_corrected: The probability after correction for multiple comparisons with the Bonferroni correction. In this correction, the p-value is multiplied by the number of comparisons performed (the number of OTUs remaining after applying the filter). 
* FDR_corrected: The probability after correction with the "false discovery rate" method. In this method, the raw p-values are ranked from low to high. Each p-value is multiplied by the number of comparisons divided by the rank. This correction is less conservative than the Bonferroni correction. The list of significant OTUs is expected to have the percent of false positives predicted by the p value.
* Category Mean Columns: Contains one column for each category reporting the mean count of the OTU in that category.
* Consensus lineage: The consensus lineage for that OTU will be listed in the last column if it was present in the input OTU table.

The correlation and longitudinal_correlation test results are output as tab 
delimited text, which can be examined in Excel. The output has the following columns:

* OTU: The name of the OTU.  
* prob: The probability that the OTU relative abundance is correlated with the category values across samples. 
* otu_values_y: a list of the values (relative abundance) of the OTU across the samples that were plotted on the y axis for the correlation.
* cat_values_x: a list of the values of the selected category that were plotted on the x axis for the correlation.
* Bonferroni_corrected: The probability after correction for multiple comparisons with the Bonferroni correction. In this correction, the p-value is multiplied by the number of comparisons performed (the number of OTUs remaining after applying the filter). 
* FDR_corrected: The probability after correction with the "false discovery rate" method. In this method, the raw p-values are ranked from low to high. Each p-value is multiplied by the number of comparisons divided by the rank. This correction is less conservative than the Bonferroni correction. The list of significant OTUs is expected to have the percent of false positives predicted by the p value.
* r: Pearson's r. This value ranges from -1 to +1, with -1 indicating a perfect negative correlation, +1 indicating a perfect positive correlation, and 0 indicating no relationship.
* Consensus lineage: The consensus lineage for that OTU will be listed in the last column if it was present in the input OTU table.

The paired_T results are output as tab delimited text that can be examined in Excel. The output has the following columns:

* OTU: The name of the OTU.
* prob: The raw probability from the paired T test
* T stat: The raw T value
* average_diff: The average difference between the before and after samples in the individuals in which the OTU was observed.
* num_pairs: The number of sample pairs (individuals) in which the OTU was observed.
* Bonferroni_corrected: The probability after correction for multiple comparisons with the Bonferroni correction. In this correction, the p-value is multiplied by the number of comparisons performed (the number of OTUs remaining after applying the filter). 
* FDR_corrected: The probability after correction with the "false discovery rate" method. In this method, the raw p-values are ranked from low to high.  Each p-value is multiplied by the number of comparisons divided by the rank.  This correction is less conservative than the Bonferroni correction. The list of significant OTUs is expected to have the percent of false positives predicted by the p value.
* Consensus lineage: The consensus lineage for that OTU will be listed in the last column if it was present in the input OTU table.




**G-test:**

Perform a G test on otu_table.biom testing OTUs for differences in the abundance across the category "Treatment":

::

	otu_category_significance.py -i otu_table.biom -m Fasting_Map.txt -s g_test -c Treatment -o single_g_test.txt

**ANOVA:**

Perform an ANOVA on otu_table.biom testing OTUs for differences in the abundance across the category "Treatment":

::

	otu_category_significance.py -i otu_table.biom -m Fasting_Map.txt -s ANOVA -c Treatment -o single_anova.txt

**ANOVA on multiple OTU tables and collate results:**

Perform an ANOVA on all OTU tables in rarefied_otu_tables testing OTUs for differences in the abundance across the category "Treatment" and collate the results into one file:

::

	otu_category_significance.py -i rarefied_otu_tables/ -m Fasting_Map.txt -s ANOVA -c Treatment -o multiple_anova.txt -w 

**ANOVA on multiple OTU tables and write out separate results files:**

Perform an ANOVA on all OTU tables in rarefied_otu_tables testing OTUs for differences in the abundance across the category "Treatment" and produce a results file for each OTU table:

::

	otu_category_significance.py -i rarefied_otu_tables/ -m Fasting_Map.txt -s ANOVA -c Treatment -o rarefied_table_results/ 


