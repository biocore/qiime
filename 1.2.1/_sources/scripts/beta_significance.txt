.. _beta_significance:

.. index:: beta_significance.py

*beta_significance.py* -- This script runs any of a set of common tests to determine if a sample is statistically significantly different from another sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The tests are conducted on each pair of samples present in the input otu table. See the unifrac tutorial online for more details (http://bmf2.colorado.edu/unifrac/tutorial.psp)


**Usage:** :file:`beta_significance.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input otu table
	-o, `-`-output_path
		Output results path
	-s, `-`-significance_test
		Significance test to use, options are 'unweighted_unifrac', 'weighted_unifrac', or 'p-test'
	
	**[OPTIONAL]**
		
	-t, `-`-tree_path
		Path to newick tree file, required for phylogenetic metrics [default: None]
	-n, `-`-num_iters
		Number of monte carlo randomizations [default: 100]


**Output:**

The script outputs a tab delimited text file with each pair of samples and a p value representing the probability that a random sample/sequence assignment will result in more dissimilar samples than the actual pair of samples.


**Example:**

Perform 100 randomizations of sample/sequence assignments, and record the probability that sample 1 is phylogenetically different from sample 2, using the unifrac monte carlo significance test. The test is run for all pairs of samples.

::

	python beta_significance.py -i otu_table.txt -t rep_set.tre -s unweighted_unifrac -o unw_sig.txt


