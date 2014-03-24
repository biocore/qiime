.. _beta_diversity:

.. index:: beta_diversity.py

*beta_diversity.py* -- Calculate beta diversity (pairwise sample dissimilarity) on one or many otu tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The input for this script is the OTU table containing the number of sequences observed in each OTU (rows) for each sample (columns). For more information pertaining to the OTU table refer to the documentation for make_otu_table. If the user would like phylogenetic beta diversity metrics using UniFrac, a phylogenetic tree must also be passed as input (see `make_phylogeny.py <./make_phylogeny.html>`_). The output of this script is a distance matrix containing a dissimilarity value for each pairwise comparison.

A number of metrics are currently supported, including unweighted and weighted UniFrac (pass the -s option to see available metrics). In general, because unifrac uses phylogenetic information, one of the unifrac metrics is recommended, as results can be vastly more useful (Hamady & Knight, 2009). Quantitative measures (e.g. weighted unifrac) are ideally suited to revealing community differences that are due to changes in relative taxon abundance (e.g., when a particular set of taxa flourish because a limiting nutrient source becomes abundant). Qualitative measures (e.g. unweighted unifrac) are most informative when communities differ primarily by what can live in them (e.g., at high temperatures), in part because abundance information can obscure significant patterns of variation in which taxa are present (Lozupone et al., 2007). Most qualitative measures are referred to here e.g. "binary_jaccard". Typically both weighted and unweighted unifrac are used.


**Usage:** :file:`beta_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input_path
		Input OTU table in biom format or input directory containing OTU tables in biom format for batch processing.
	-r, `-`-rows
		Compute for only these rows of the distance matrix. User should pass a list of sample names (e.g. "s1,s3") [default: None; full n x n matrix is generated]
	-o, `-`-output_dir
		Output directory. One will be created if it doesn't exist.
	-m, `-`-metrics
		Beta-diversity metric(s) to use. A comma-separated list should be provided when multiple metrics are specified. [default: unweighted_unifrac,weighted_unifrac]
	-s, `-`-show_metrics
		Show the available beta-diversity metrics and exit. Metrics starting with "binary..." specifies that a metric is qualitative, and considers only the presence or absence of each taxon [default: False]
	-t, `-`-tree_path
		Input newick tree filepath, which is required when phylogenetic metrics are specified. [default: None]
	-f, `-`-full_tree
		By default, tips not corresponding to OTUs in the OTU table are removed from the tree for diversity calculations. Pass to skip this step if you're already passing a minimal tree. Beware with "full_tree" metrics, as extra tips in the tree change the result


**Output:**

Each file in the input directory should be an otu table, and the output of `beta_diversity.py <./beta_diversity.html>`_ is a folder containing text files, each a distance matrix between samples corresponding to an input otu table.


**Single File Beta Diversity (non-phylogenetic):**

To perform beta diversity (using e.g. euclidean distance) on a single OTU table, where the results are output to beta_div/, use the following command:

::

	beta_diversity.py -i otu_table.biom -m euclidean -o beta_div

**Single File Beta Diversity (phylogenetic):**

In the case that you would like to perform beta diversity using a phylogenetic metric (e.g. weighted_unifrac), you can use the following command:

::

	beta_diversity.py -i otu_table.biom -m weighted_unifrac -o beta_div/ -t rep_set.tre

**Multiple File (batch) Beta Diversity (phylogenetic):**

To perform beta diversity on multiple OTU tables (e.g., resulting files from `multiple_rarefactions.py <./multiple_rarefactions.html>`_), specify an input directory (e.g. otu_tables/) as shown by the following command:

::

	beta_diversity.py -i otu_tables/ -m weighted_unifrac -o beta_div/ -t rep_set.tre


