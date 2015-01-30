.. _conditional_uncovered_probability:

.. index:: conditional_uncovered_probability.py

*conditional_uncovered_probability.py* -- Calculate the conditional uncovered probability on each sample in an otu table.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script calculates the conditional uncovered probability for each sample
in an OTU table. It uses the methods introduced in Lladser, Gouet, and Reeder,
"Extrapolation of Urn Models via Poissonization: Accurate Measurements of the
Microbial Unknown" PLoS 2011.

Specifically, it computes a point estimate and a confidence interval using two
different methods. Thus it can happen that the PE is actually outside of the
CI.

We only provide the ability to generate 95% (alpha=0.95) CIs. The CIs are ULCL
CIs; they provide an upper and lower bound, where the lower bound is
conservative. The CIs are constructed using an upper-to-lower bound ratio of
10.

The CI method requires precomputed constants that depend on the lookahead. We
only provide constants for r=3..25,30,40,50.




**Usage:** :file:`conditional_uncovered_probability.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input_path
		Input OTU table filepath. [default: None]
	-o, `-`-output_path
		Output filepath to store the predictions. [default: None]
	-r, `-`-look_ahead
		Number of unobserved, new colors necessary for prediction. [default: 25]
	-m, `-`-metrics
		CUP metric(s) to use. A comma-separated list should be provided when multiple metrics are specified. [default: lladser_pe,lladser_ci]
	-s, `-`-show_metrics
		Show the available CUP metrics and exit.


**Output:**

The resulting file(s) is a tab-delimited text file, where the columns
correspond to estimates of the cond. uncovered probability and the rows
correspond to samples. The output file is compatible with the alpha_diversity
output files and thus could be tied into the rarefaction workflow.

Example Output:

====== ======= ============= ================
\      PE      Lower Bound   Upper Bound
====== ======= ============= ================
PC.354 0.111   0.0245        0.245
PC.124 0.001   0.000564      0.00564
====== ======= ============= ================




**Default case:**

To calculate the cond. uncovered probability with the default values, you can use the following command:

::

	conditional_uncovered_probability.py -i otu_table.biom -o cup.txt

**Change lookahead:**

To change the accuracy of the prediction change the lookahead value. Larger values of r lead to more precise predictions, but might be unfeasable for small samples. For deeply sequenced samples, try increasing r to 50:

::

	conditional_uncovered_probability.py -i otu_table.biom -o cup_r50.txt -r 50


