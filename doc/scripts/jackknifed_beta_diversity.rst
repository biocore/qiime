.. _jackknifed_beta_diversity:

.. index:: jackknifed_beta_diversity.py

*jackknifed_beta_diversity.py* -- A workflow script for performing jackknifed UPGMA clustering and building jackknifed Emperor PCoA plots.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To directly measure the robustness of individual UPGMA clusters and clusters in PCoA plots, one can perform jackknifing (repeatedly resampling a subset of the available data from each sample).


**Usage:** :file:`jackknifed_beta_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input OTU table in biom format [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	-e, `-`-seqs_per_sample
		Number of sequences to include in each jackknifed subset [REQUIRED]
	-m, `-`-mapping_fp
		Path to the mapping file [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-tree_fp
		Path to the tree file [default: None; REQUIRED for phylogenetic measures]
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	`-`-master_tree
		Method for computing master trees in jackknife analysis. "consensus": consensus of trees from jackknifed otu tables.  "full": tree generated from input (unsubsambled) otu table.  [default: consensus]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]


**Output:**

This scripts results in several distance matrices (from `beta_diversity.py <./beta_diversity.html>`_), several rarified OTU tables (from `multiple_rarefactions_even_depth.py <./multiple_rarefactions_even_depth.html>`_), several UPGMA trees (from `upgma_cluster.py <./upgma_cluster.html>`_), a supporting file and newick tree with support values (from `tree_compare.py <./tree_compare.html>`_), and Emperor PCoA plots.


**Example:**

These steps are performed by the following command: Compute beta diversity distance matrix from otu table (and tree, if applicable); build rarefied OTU tables by evenly sampling to the specified depth (-e); build UPGMA tree from full distance matrix; compute distance matrics for rarefied OTU tables; build UPGMA trees from rarefied OTU table distance matrices; build a consensus tree from the rarefied UPGMA trees; compare rarefied OTU table distance matrix UPGMA trees to either (full or consensus) tree for jackknife support of tree nodes; perform principal coordinates analysis on distance matrices generated from rarefied OTU tables; generate Emperor PCoA plots with jackknifed support.



::

	jackknifed_beta_diversity.py -i otu_table.biom -o bdiv_jk100 -e 100 -m Fasting_Map.txt -t rep_set.tre


