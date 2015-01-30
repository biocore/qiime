.. _transform_coordinate_matrices:

.. index:: transform_coordinate_matrices.py

*transform_coordinate_matrices.py* -- Transform two or more coordinate matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script transforms two or more coordinate matrices (e.g., the output of `principal_coordinates.py <./principal_coordinates.html>`_) using procrustes analysis to minimize the distances between corresponding points. The first coordinate matrix provided is treated as the reference, and all other coordinate matrices are transformed to minimize distances to the reference points. Monte Carlo simulations can additionally be performed (-r random trials are run) to estimate the probability of seeing an M^2 value as extreme as the actual M^2.


**Usage:** :file:`transform_coordinate_matrices.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		Comma-separated list of input coordinate matrices
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-r, `-`-random_trials
		Number of random permutations of matrix2 to perform. [default: (no Monte Carlo analysis performed)]
	-d, `-`-num_dimensions
		Number of dimensions to include in output matrices [default: 3]
	-s, `-`-sample_id_map_fps
		If sample id maps are provided, there must be exactly one fewer files here than there are coordinate matrices (as each nth sample id map will provide the mapping from the first input coordinate matrix to the n+1th coordinate matrix) [default: None]
	`-`-store_trial_details
		Store PC matrices for individual trials [default: False]


**Output:**

Two transformed coordinate matrices corresponding to the two input coordinate matrices, and (if -r was specified) a text file summarizing the results of the Monte Carlo simulations.


**Write the transformed procrustes matrices to file:**

::

	transform_coordinate_matrices.py -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt -o procrustes_output

**Generate transformed procrustes matrices and monte carlo p-values for two principal coordinate matrices:**

::

	transform_coordinate_matrices.py -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt -o mc_procrustes_output_2 -r 1000

**Generate transformed procrustes matrices and monte carlo p-values for four principal coordinate matrices:**

::

	transform_coordinate_matrices.py -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt,euclidean_pc.txt,bray_curtis_pc.txt -o mc_procrustes_output_4 -r 1000

**Generate transformed procrustes matrices and monte carlo p-values for three principal coordinate matrices where the sample ids must be mapped between matrices:**

::

	transform_coordinate_matrices.py -i s1_pc.txt,s2_pc.txt,s3_pc.txt -s s1_s2_map.txt,s1_s3_map.txt -o mc_procrustes_output_3 -r 1000


