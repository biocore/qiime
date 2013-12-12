.. _identify_paired_differences:

.. index:: identify_paired_differences.py

*identify_paired_differences.py* -- Generate plots and stats to test for change in some data point(s) with a state change on a per-individual basis.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script provides a framework for paired-difference testing (i.e., analysis of data generated under a pre/post experimental design). In a pre/post experimental design, individuals are sampled before and after some 'treatment'. This code plots differences in values in the sample metadata (i.e., the mapping file) or observation counts in a BIOM table, and runs a (Bonferroni-corrected) one sample t-test on each sample metadata category or BIOM observation to determine if the mean of each distribution of pre/post differences differs from zero. If 'None' appears for the t score and p-values, this often means that the distribution of differences contained no variance, so the t-test could not be run. This can happen, for example, if the value passed for --valid_states is so restrictive that only a single sample is retained for analysis.


**Usage:** :file:`identify_paired_differences.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		The input metadata map filepath
	-o, `-`-output_dir
		Directory where output files should be saved
	-t, `-`-state_category
		The mapping file column name to plot change over (usually has values like "pre-treatment" and "post-treatment")
	-x, `-`-state_values
		Ordered list of state values to test change over (defines direction of graphs, generally something like "pre-treatment,post-treatment"). currently limited to two states.
	-c, `-`-individual_id_category
		The mapping file column name containing each individual's identifier (usually something like "personal_identifier")
	
	**[OPTIONAL]**
		
	`-`-ymin
		Set the minimum y-value across plots [default: determined on a per-plot basis]
	`-`-ymax
		Set the maximum y-value across plots [default: determined on a per-plot basis]
	`-`-metadata_categories
		Ordered list of the mapping file column names to test for paired differences (usually something like "StreptococcusAbundance,Phylogenetic Diversity") [default: None]
	`-`-observation_ids
		Ordered list of the observation ids to test for paired differences if a biom table is provided (usually something like "otu1,otu2") [default: compute paired differences for all observation ids]
	-b, `-`-biom_table_fp
		Path to biom table to use for computing paired differences [default: None]
	-s, `-`-valid_states
		String describing samples that should be included based on their metadata (e.g. 'TreatmentResponse:Improved') [default: all samples are included in analysis]
	`-`-line_color
		Color of lines in plots, useful if generating multiple plots in different runs of this script to overlay on top of one another. these can be specified as matplotlib color names, or as html hex strings [default: black]


**Output:**

The output of this script is plots of pre/post differences and associated statistics.


**Generate plots and stats for one category from the mapping file where the y-axis should be consistent across plots and the lines in the plots should be light blue.:**

::

	identify_paired_differences.py -m map.txt --metadata_categories 'Streptococcus Abundance' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_results --ymin 0 --ymax 60 --line_color '#eeefff'

**Generate plots and stats for three categories from the mapping file.:**

::

	identify_paired_differences.py -m map.txt --metadata_categories 'Streptococcus Abundance,Phylogenetic Diversity,Observed OTUs' --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o taxa_and_alpha_results

**Generate plots for all observations in a biom file:**

::

	identify_paired_differences.py -m map.txt -b otu_table.biom --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o otu_results

**Generate plots for all observations in a biom file, but only including samples from individuals whose 'TreatmentResponse' was 'Improved' (as defined in the mapping file).:**

::

	identify_paired_differences.py -m map.txt -b otu_table.biom --state_category TreatmentState --state_values Pre,Post --individual_id_category PersonalID -o otu_results_improved_only --valid_states TreatmentResponse:Improved


