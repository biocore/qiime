.. _alpha_rarefaction:

.. index:: alpha_rarefaction.py

*alpha_rarefaction.py* -- A workflow script for performing alpha rarefaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


The steps performed by this script are:

1. Generate rarefied OTU tables;

2. Compute alpha diversity metrics for each rarefied OTU table;

3. Collate alpha diversity results;

4. Generate alpha rarefaction plots.



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
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]


**Output:**

The results of this script is a folder ("rare1/") containing rarefied otu tables, alpha-diversity for each otu table, a file containing the results from collating the alpha-diversity results and a folder containing the rarefaction plots.


**Example:**

::

	alpha_rarefaction.py -o rare1 -i otu_table.txt -t inseqs1_rep_set.tre -m inseqs1_mapping.txt -p custom_parameters.txt


