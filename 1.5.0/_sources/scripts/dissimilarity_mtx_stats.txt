.. _dissimilarity_mtx_stats:

.. index:: dissimilarity_mtx_stats.py

*dissimilarity_mtx_stats.py* -- Calculate mean, median and standard deviation from a set of distance matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script reads in all (dis)similarity matrices from an input directory (input_dir), then calculates and writes the mean, median, standdard deviation (stdev) to an output folder.

The input_dir must contain only (dis)similarity matrices, and only those you wish to perform statistical analyses on.


**Usage:** :file:`dissimilarity_mtx_stats.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Path to input directory
	-o, `-`-output_dir
		Path to store result files


**Output:**

The outputs are in distance matrix format, where each value is the mean, median, or stdev of that element in all the input distance matrices


**Example:**

This examples takes the "dists/" directory as input and returns the results in the "dist_stats/" directory:

::

	dissimilarity_mtx_stats.py -i dists/ -o dist_stats/


