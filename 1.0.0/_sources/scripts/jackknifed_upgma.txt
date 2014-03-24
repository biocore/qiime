.. _jackknifed_upgma:

.. index:: jackknifed_upgma.py

*jackknifed_upgma.py* -- A workflow script for performing jackknifed UPGMA clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To directly measure the robustness of individual UPGMA clusters, one can perform jackknifing (repeatedly resampling a subset of the available data from each sample). This process will utilize the following script: `beta_diversity.py <./beta_diversity.html>`_, `multiple_rarefactions.py <./multiple_rarefactions.html>`_, `upgma_cluster.py <./upgma_cluster.html>`_ and `tree_compare.py <./tree_compare.html>`_).


**Usage:** :file:`jackknifed_upgma.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input fasta file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	-p, `-`-parameter_fp
		Path to the parameter file [REQUIRED]
	-e, `-`-seqs_per_sample
		Number of sequences to include in each jackknifed subset [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-tree_fp
		Path to the tree file [default: None; REQUIRED for phylogenetic measures]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]


**Output:**

This scripts results in several distance matrices (from `beta_diversity.py <./beta_diversity.html>`_), several rarified otu tables (from `multiple_rarefactions.py <./multiple_rarefactions.html>`_) several UPGMA trees (from `upgma_cluster.py <./upgma_cluster.html>`_) and a supporting file and newick tree with support values (from `tree_compare.py <./tree_compare.html>`_).


**Example:**

These steps are performed by the following command:

1. Compute beta diversity distance matrix from otu table (and tree, if applicable)

2. Build rarefied OTU tables;

3. Build UPGMA tree from full distance matrix;

4. Compute distance matrics for rarefied OTU tables; 

5. Build UPGMA trees from rarefied OTU table distance matrices;

6. Compare rarefied OTU table distance matrix UPGMA trees to tree full UPGMA tree and write support file and newick tree with support values as node labels.



::

	jackknifed_upgma.py -i inseqs1_otu_table.txt -t inseqs1_rep_set.tre -p custom_parameters_jack.txt -o wf_jack -e 5 -v


