.. _conditional_uncovered_probability:

.. index:: conditional_uncovered_probability.py

*conditional_uncovered_probability.py* -- Calculate the conditional uncovered probability on each sample in an otu table.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script calculates the conditional uncovered probability for each sample in an OTU table. It uses the methods introduced in Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization: Accurate Measurements of the Microbial Unknown" PLoS 2011. 

Specifically, it computes a point estimate and a confidence interval using two different methods. Thus it can happen that the PE is actually outside of the CI. 

The CI method requires precomputed constants that depend on the lookahead, the upper-to-lower bound ratio and the desired confidence.
We only provide these constants for some frequently used combinations. These are (alpha:0.95, r=1..25)) for the the L and U interval types, and (alpha:0.9, 0.95, 0.99; f=10;  r=3..25,30,40,50). Also, there are a few hand picked special cases:

f=2 and r=50 and alpha=0.95
f=2 and r=33 and alpha=0.95
f=1.5 and r=100 and alpha=0.95
f=1.5 and r=94 and alpha=0.95
f=2.5 and r=19 and alpha=0.95




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
	-c, `-`-ci_type
		Type of confidence interval.  Choice of ULCL, ULCU, U, L [default: ULCL]
	-a, `-`-alpha
		Desired confidence level for CI prediction. [default: 0.95]
	-f, `-`-f_ratio
		Upper to lower bound ratio for CI prediction. [default: 10.0]
	-m, `-`-metrics
		CUP metric(s) to use. A comma-separated list should be provided when multiple metrics are specified. [default: lladser_pe,lladser_ci]
	-s, `-`-show_metrics
		Show the available CUP metrics and exit.


**Output:**

The resulting file(s) is a tab-delimited text file, where the columns correspond to estimates of the cond. uncovered probability and the rows correspond to samples. The output file is compatible with the alpha_diversity output files and thus could be tied into thes rarefaction workflow.

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

**Change the interval type:**

To change the confidence interval type to a lower bound prediction, while the upper bound is set to 1 use: 

::

	conditional_uncovered_probability.py -i otu_table.biom -o cup_lower_bound.txt -c L


