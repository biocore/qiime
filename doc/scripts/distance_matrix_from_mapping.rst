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
		String containing the name of the column in the mapping file, e.g. 'DOB'. If you pass two colums separated by a comma (e.g. 'Latitude,Longitud') the script will calculate the Vincenty formula (WGS-84) for distance between two Latitude/Longitude points.
	
	**[OPTIONAL]**
		
	-o, `-`-output_fp
		Output directory. One will be created if it doesn't exist. [default=map_distance_matrix.txt]


**Output:**

The output of `distance_matrix_from_mapping.py <./distance_matrix_from_mapping.html>`_ is a file containing a distance matrix between rows corresponding to a pair of columns in a mapping file.


**Pairwise dissimilarity:**

To calculate the distance matrix (using euclidean distance) on a column of the mapping file, where the results are output to DOB.txt, use the following command:

::

	distance_matrix_from_mapping.py -i Fasting_Map.txt -c DOB

**Pairwise dissimilarity using the Vincenty formula for distance between two Latitude/Longitude points:**

To calculate the distance matrix (using Vincenty formula) on a column of the mapping file, where the results are output to lat_long.txt, use the following command:

::

	distance_matrix_from_mapping.py -i lat_long.txt -c Latitute,Longitude -o lat_long_dtx_matrix.txt


