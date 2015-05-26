.. _alpha_diversity:

.. index:: alpha_diversity.py

*alpha_diversity.py* -- Calculate alpha diversity on each sample in an otu table, using a variety of alpha diversity metrics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script calculates alpha diversity, or within-sample diversity, using an otu table. The QIIME pipeline allows users to conveniently calculate more than two dozen different diversity metrics. The full list of available metrics is available by passing the option -s to the script `alpha_diversity.py <./alpha_diversity.html>`_. Every metric has different strengths and limitations - technical discussion of each metric is readily available online and in ecology textbooks, but is beyond the scope of this document.


**Usage:** :file:`alpha_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input_path
		Input OTU table filepath or input directory containing OTU tables for batch processing. [default: None]
	-o, `-`-output_path
		Output distance matrix filepath or output directory to store distance matrices when batch processing. [default: None]
	-m, `-`-metrics
		Alpha-diversity metric(s) to use. A comma-separated list should be provided when multiple metrics are specified. [default: PD_whole_tree,chao1,observed_species]
	-s, `-`-show_metrics
		Show the available alpha-diversity metrics and exit.
	-t, `-`-tree_path
		Input newick tree filepath. [default: None; REQUIRED for phylogenetic metrics]


**Output:**

The resulting file(s) is a tab-delimited text file, where the columns correspond to alpha diversity metrics and the rows correspond to samples and their calculated diversity measurements. When a folder is given as input (-i), the script processes every otu table file in the given folder, and creates a corresponding file in the output directory.

Example Output:

====== ======= ============= ================
\      simpson PD_whole_tree observed_species
====== ======= ============= ================
PC.354 0.925   2.83739       16.0
PC.355 0.915   3.06609       14.0
PC.356 0.945   3.10489       19.0
PC.481 0.945   3.65695       19.0
PC.593 0.91    3.3776        15.0
PC.607 0.92    4.13397       16.0
PC.634 0.9     3.71369       14.0
PC.635 0.94    4.20239       18.0
PC.636 0.925   3.78882       16.0
====== ======= ============= ================



**Single File Alpha Diversity Example (non-phylogenetic):**

To perform alpha diversity (e.g. chao1) on a single OTU table, where the results are output to "alpha_div.txt", you can use the following command:

::

	alpha_diversity.py -i otu_table.biom -m chao1 -o adiv_chao1.txt

**Single File Alpha Diversity Example (phylogenetic):**

In the case that you would like to perform alpha diversity using a phylogenetic metric (e.g. PD_whole_tree), you can use the following command:

::

	alpha_diversity.py -i otu_table.biom -m PD_whole_tree -o adiv_pd.txt -t rep_set.tre

**Single File Alpha Diversity Example with multiple metrics:**

You can use the following idiom to run multiple metrics at once (comma-separated):

::

	alpha_diversity.py -i otu_table.biom -m chao1,PD_whole_tree -o adiv_chao1_pd.txt -t rep_set.tre

**Multiple File (batch) Alpha Diversity:**

To perform alpha diversity on multiple OTU tables (e.g.: rarefied otu tables resulting from `multiple_rarefactions.py <./multiple_rarefactions.html>`_), specify an input directory instead of a single otu talbe, and an output directory (e.g. "alpha_div_chao1_PD/") as shown by the following command:

::

	alpha_diversity.py -i otu_tables/ -m chao1,PD_whole_tree -o adiv_chao1_pd/ -t rep_set.tre


