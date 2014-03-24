.. _transform_coordinate_matrices:

.. index:: transform_coordinate_matrices.py

*transform_coordinate_matrices.py* -- Transform 2 coordinate matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script transformms 2 coordinate matrices (e.g., the output of `principal_coordinates.py <./principal_coordinates.html>`_) using procrustes analysis to minimize the distances between corresponding points. Monte Carlo simulations can additionally be performed (-r random trials are run) to estimate the probability of seeing an M^2 value as extreme as the actual M^2.


**Usage:** :file:`transform_coordinate_matrices.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		Comma-separated input files
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-r, `-`-random_trials
		Number of random permutations of matrix2 to perform.  [default: (no Monte Carlo analysis performed)]
	-d, `-`-num_dimensions
		Number of dimensions to include in output matrices [default: Consider all dimensions]
	-s, `-`-sample_id_map_fp
		Map of original sample ids to new sample ids [default: None]


**Output:**

Two transformed coordinate matrices corresponding to the two input coordinate matrices, and (if -r was specified) a text file summarizing the results of the Monte Carlo simulations.


**Generate monte carlo p-values:**

::

	transform_coordinate_matrices.py -r 1000 -i weighted_unifrac_coords.txt,unweighted_unifrac_coords.txt

**Write the transformed procrustes matrices to file:**

::

	transform_coordinate_matrices.py -o out/ -i weighted_unifrac_coords.txt,unweighted_unifrac_coords.txt


