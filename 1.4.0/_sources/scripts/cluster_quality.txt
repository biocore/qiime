.. _cluster_quality:

.. index:: cluster_quality.py

*cluster_quality.py* -- compute the quality of a cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The input is a distance matrix (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_).


**Usage:** :file:`cluster_quality.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input distance matrix file
	-m, `-`-map
		Mapping file
	-c, `-`-category
		Column of mapping file delimiting clusters
	
	**[OPTIONAL]**
		
	-o, `-`-output_path
		Output path, prints to stdout if omitted
	-s, `-`-short
		Print only the ratio of mean dissimilarities between/within clusters instead of more detailed output
	`-`-metric
		Choice of quality metric to apply. Currently only one option exists, the ratio of mean(distances between samples from different clusters) to mean(distances between samples from the same cluster) Default: ratio


**Output:**

The output is either a single number (with -s), or a more detailed output of the similarity between and within clusters.


**cluster quality based on the treatment category:**

to compute the quality of clusters, and print to stdout, use the following idiom:

::

	cluster_quality.py -i unweighted_unifrac_distance_matrix.txt -m Fasting_Map.txt -c Treatment


