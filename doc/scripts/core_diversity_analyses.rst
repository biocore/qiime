.. _core_diversity_analyses:

.. index:: core_diversity_analyses.py

*core_diversity_analyses.py* -- A workflow for running a core set of QIIME diversity analyses.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script plugs several QIIME diversity analyses together to form a basic workflow beginning with a BIOM table, mapping file, and optional phylogenetic tree. 

The included scripts are those run by the workflow scripts `alpha_rarefaction.py <./alpha_rarefaction.html>`_, `beta_diversity_through_plots.py <./beta_diversity_through_plots.html>`_, `summarize_taxa_through_plots.py <./summarize_taxa_through_plots.html>`_, plus the (non-workflow) scripts `make_distance_boxplots.py <./make_distance_boxplots.html>`_, `compare_alpha_diversity.py <./compare_alpha_diversity.html>`_, and `otu_category_significance.py <./otu_category_significance.html>`_. To update parameters to the workflow scripts, you should pass the same parameters file that you would pass if calling the workflow script directly.



**Usage:** :file:`core_diversity_analyses.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_biom_fp
		The input biom file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	-m, `-`-mapping_fp
		The mapping filepath [REQUIRED]
	-e, `-`-sampling_depth
		Sequencing depth to use for even sub-sampling and maximum rarefaction depth. You should review the output of `biom summarize-table <http://biom-format.org/documentation/summarizing_biom_tables.html>`_ to decide on this value.
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. For more information, see www.qiime.org/documentation/qiime_parameters_files.html [if omitted, default values will be used]
	-a, `-`-parallel
		Run in parallel where available. Specify number of jobs to start with -O or in the parameters file. [default: False]
	`-`-nonphylogenetic_diversity
		Apply non-phylogenetic alpha (chao1 and observed_species) and beta (bray_curtis) diversity calculations. This is useful if, for example, you are working with non-amplicon BIOM tables, or if a reliable tree is not available (e.g., if you're  working with ITS amplicons) [default: False]
	`-`-suppress_taxa_summary
		Suppress generation of taxa summary plots. [default: False]
	`-`-suppress_beta_diversity
		Suppress beta diversity analyses. [default: False]
	`-`-suppress_alpha_diversity
		Suppress alpha diversity analyses. [default: False]
	`-`-suppress_otu_category_significance
		Suppress OTU/category significance analysis. [default: False]
	-t, `-`-tree_fp
		Path to the tree file if one should be used. [default: no tree will be used]
	-c, `-`-categories
		The metadata category or categories to compare (i.e., column headers in the mapping file) for categorical analyses. These should be passed  as a comma-separated list. [default: None; do not perform categorical analyses]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging or recovering from failed runs. [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 2]


**Output:**




Run diversity analyses at 20 sequences/sample, with categorical analyses focusing on the SampleType and day categories. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	core_diversity_analyses.py -i $PWD/otu_table.biom -o $PWD/core_output -m $PWD/map.txt -c SampleType,day -t $PWD/rep_set.tre -e 20


