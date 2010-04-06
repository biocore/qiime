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
		Input path.  directory for batch processing, filename for single file operation
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation
	-m, `-`-metrics
		Metrics to use, comma delimited
	-s, `-`-show_metrics
		Show available alpha diversity metrics and quit
	-t, `-`-tree_path
		Path to newick tree file, required for phylogenetic metrics [default: None]


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



**Single File Alpha Diversity Example:**

To perform alpha diversity (e.g. chao1) on a single OTU table, where the results are output to "alpha_div.txt", you can use the following command:

::

	alpha_diversity.py -i otu_table.txt -m chao1 -o alpha_div.txt

Note: Since this is a non-phylogenetic metric, the tree does not need to be supplied.

In the case that you would like to perform alpha diversity using a phylogenetic metric (e.g. PD_whole_tree), you can use the following command:

::

	alpha_diversity.py -i otu_table.txt -m PD_whole_tree -o alpha_div.txt -t repr_set.tre

You can use the following idiom to run multiple metrics at once (comma-separated):

::

	alpha_diversity.py -i otu_table.txt -m chao1,PD_whole_tree -o alpha_div.txt -t repr_set.tre

**Multiple File (batch) Alpha Diversity:**

To perform alpha diversity on multiple OTU tables (e.g.: rarefied otu tables resulting from `multiple_rarefactions.py <./multiple_rarefactions.html>`_), specify an input directory instead of a single otu talbe, and an output directory (e.g. "alpha_div_chao1_PD/") as shown by the following command:

::

	alpha_diversity.py -i rarefaction_tables/ -m chao1,PD_whole_tree -o alpha_div_chao1_PD/ -t repr_set.tre


