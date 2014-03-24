.. _beta_diversity:

.. index:: beta_diversity.py

*beta_diversity.py* -- Calculate beta diversity (pairwise sample dissimilarity) on one or many otu tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The input for this script is the OTU table containing the number of sequences observed in each OTU (rows) for each sample (columns). For more information pertaining to the OTU table refer to the documentation for make_otu_table. If the user would like phylogenetic beta diversity metrics using UniFrac, a phylogenetic tree must also be passed as input (see `make_phylogeny.py <./make_phylogeny.html>`_). The output of this script is a distance matrix containing a dissimilarity value for each pairwise comparison.

A number of metrics are currently supported, including unweighted and weighted UniFrac (pass the -s option to see available metrics). In general, because unifrac uses phylogenetic information, one of the unifrac metrics is recommended, as results can be vastly more useful (Hamady & Knight, 2009). Weighted unifrac is sensitive to factors affecting the relative abundance of different taxa, such as salinity or pH while unweighted unifrac considers only the presence or absence of taxa, and is sensitive to factors such as nutrient availability (Costello, personal communication). Typically both weighted and unweighted unifrac are applied.


**Usage:** :file:`beta_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input_path
		Input path: otu table, or dir of otu tables for batch mode
	-o, `-`-output_dir
		Output directory, will be created if doesn't exist
	-m, `-`-metrics
		Metrics to use, comma delimited if >1 metric, no spaces
	-s, `-`-show_metrics
		Show available beta diversity metrics and quit
	-t, `-`-tree_path
		Path to newick tree file, required for phylogenetic metrics [default: None]


**Output:**

Each file in the input directory should be an otu table, and the output of `beta_diversity.py <./beta_diversity.html>`_ is a folder containing text files, each a distance matrix between samples corresponding to an input otu table.


**Single File Beta Diversity:**

To perform beta diversity (using e.g. euclidean distance) on a single OTU table, where the results are output to beta_div.txt, use the following command:

::

	beta_diversity.py -i otu_table.txt -m euclidean -o beta_div/

Note: Since this is a non-phylogenetic metric, the tree does not need to be supplied.

In the case that you would like to perform beta diversity using a phylogenetic metric (e.g. weighted_unifrac), you can use the following command:

::

	beta_diversity.py -i otu_table.txt -m weighted_unifrac -o beta_div/ -t repr_set.tre

**Multiple File (batch) Beta Diversity:**

To perform beta diversity on multiple OTU tables (resulting files from `multiple_rarefactions.py <./multiple_rarefactions.html>`_), specify an input directory (e.g. rarefaction_tables/) as shown by the following command:

::

	beta_diversity.py -i rarefaction_tables/ -m weighted_unifrac -o beta_div/ -t repr_set.tre


