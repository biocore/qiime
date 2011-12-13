.. _beta_diversity_through_plots:

.. index:: beta_diversity_through_plots.py

*beta_diversity_through_plots.py* -- A workflow script for computing beta diversity distance matrices and the corresponding 3D plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script will perform beta diversity, principal coordinate anlalysis, and generate a preferences file along with 3D PCoA Plots.



**Usage:** :file:`beta_diversity_through_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input fasta file [REQUIRED]
	-m, `-`-mapping_fp
		Path to the mapping file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-tree_fp
		Path to the tree file [default: None; REQUIRED for phylogenetic measures]
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	`-`-color_by_all_fields
		Plots will have coloring for all mapping fields [default: False; only include fields with greater than one value and fewer values than the number of samples]
	-c, `-`-histogram_categories
		Mapping fields to use when plotting distance histograms [default: None]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-e, `-`-seqs_per_sample
		Depth of coverage for even sampling [default: None]
	`-`-suppress_2d_plots
		Do not generate 2D plots [default: False]
	`-`-suppress_3d_plots
		Do not generate 3D plots [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]


**Output:**

This script results in a distance matrix (from `beta_diversity.py <./beta_diversity.html>`_), a principal coordinates file (from `principal_coordinates.py <./principal_coordinates.html>`_), a preferences file (from `make_prefs_file.py <./make_prefs_file.html>`_) and  folder containing the resulting 3d PCoA plots (as an html from `make_3d_plots.py <./make_3d_plots.html>`_).


**Example:**

The following steps are performed by the command below:

1. Compute a beta diversity distance matrix;

2. Peform a principle coordinates analysis on the result of Step 1;

3. Generate a 3D prefs file for optimized coloring of continuous variables;

4. Generate a 3D plot for all mapping fields with colors optimized for continuous data;

5. Generate a 3D plot for all mapping fields with colors optimized for discrete data.


::

	beta_diversity_through_3d_plots.py -i otu_table.txt -o bdiv1 -t inseqs1_rep_set.tre -m inseqs1_mapping.txt -p custom_parameters.txt


