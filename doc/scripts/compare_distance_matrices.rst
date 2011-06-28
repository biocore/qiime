.. _compare_distance_matrices:

.. index:: compare_distance_matrices.py

*compare_distance_matrices.py* -- Script for computing Mantel correlations between as set of distance matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`compare_distance_matrices.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dms
		The input distance matrices, comma-separated
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-n, `-`-num_iterations
		The number of iterations to perform
	-s, `-`-sample_id_map_fp
		Map of original sample ids to new sample ids [default: None]


**Output:**




Perform Mantel test on all pairs of four distance matrices, including 1000 Monte Carlo iterations. Write the output to mantel_out.txt.

::

	mantel.py -i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt,weighted_unifrac_even100_dm.txt,unweighted_unifrac_even100_dm.txt -o mantel_out.txt -n 1000


