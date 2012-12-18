.. _distance_matrix_from_mapping:

.. index:: distance_matrix_from_mapping.py

*distance_matrix_from_mapping.py* -- Calculate the pairwise dissimilarity on one column of a mappping file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The input for this script is a mapping file and the name of a column, it has to be numeric, from which a distance matrix will be created. The output of this script is a distance matrix containing a dissimilarity value for each pairwise comparison.

As this is a univariate procedure only one metric is supported: d = c-b.


**Usage:** :file:`distance_matrix_from_mapping.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Mapping filepath.
	-c, `-`-column
		String containing the name of the column in the mapping file, e.g. 'DOB'
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Output directory. One will be created if it doesn't exist. [default=map_distance_matrix]


**Output:**

The output of `distance_matrix_from_mapping.py <./distance_matrix_from_mapping.html>`_ is a file containing a distance matrix between rows corresponding to a column in a mapping file.


**Pairwise dissimilarity:**

To calculate the distance matrix (using euclidean distance) on a column of the mapping file, where the results are output to DOB.txt, use the following command:

::

	distance_matrix_from_mapping.py -i Fasting_Map.txt -c DOB


