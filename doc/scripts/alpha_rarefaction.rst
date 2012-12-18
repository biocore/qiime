.. _alpha_rarefaction:

.. index:: alpha_rarefaction.py

*alpha_rarefaction.py* -- A workflow script for performing alpha rarefaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


The steps performed by this script are: Generate rarefied OTU tables; compute alpha diversity metrics for each rarefied OTU table; collate alpha diversity results; and generate alpha rarefaction plots.


**Usage:** :file:`alpha_rarefaction.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input otu table [REQUIRED]
	-m, `-`-mapping_fp
		Path to the mapping file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	-n, `-`-num_steps
		Number of steps (or rarefied OTU table sizes) to make between min and max counts [default: 10]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-t, `-`-tree_fp
		Path to the tree file [default: None; REQUIRED for phylogenetic measures]
	`-`-min_rare_depth
		The lower limit of rarefaction depths [default: 10]
	-e, `-`-max_rare_depth
		The upper limit of rarefaction depths [default: median sequence/sample count]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 4]


**Output:**

The primary interface for the results will be OUTPUT_DIR/alpha_rarefaction_plots/rarefaction_plots.html where OUTPUT_DIR is the value you specify with -o.  You can open this in a web browser for interactive alpha rarefaction plots.


**Example:**

Given an OTU table, a phylogenetic tree, a mapping file, and a max sample depth, compute alpha rarefaction plots for the PD, observed species and chao1 metrics. To specify alternative metrics pass a parameter file via -p. We generally recommend that the max depth specified here (-e) is the same as the even sampling depth provided to beta_diversity_through_plots (also -e). 

::

	alpha_rarefaction.py -i otu_table.biom -o arare_max100/ -t rep_set.tre -m Fasting_Map.txt -e 100


