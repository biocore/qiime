.. _compare_distance_matrices:

.. index:: compare_distance_matrices.py

*compare_distance_matrices.py* -- Computes Mantel correlation tests between sets of distance matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script compares two or more distance/dissimilarity matrices for correlation by providing the Mantel, partial Mantel, and Mantel correlogram matrix correlation tests.

The Mantel test will test the correlation between two matrices. The data often represents the "distance" between objects or samples.

The partial Mantel test is a first-order correlation analysis that utilizes three distance (dissimilarity) matrices. This test builds on the traditional Mantel test which is a procedure that tests the hypothesis that distances between the objects within a given matrix are linearly independent of the distances withing those same objects in a separate matrix. It builds on the traditional Mantel test by adding a third "control" matrix.

Mantel correlogram produces a plot of distance classes versus Mantel statistics. Briefly, an ecological distance matrix (e.g. UniFrac distance matrix) and a second distance matrix (e.g. spatial distances, pH distances, etc.) are provided. The second distance matrix has its distances split into a number of distance classes (the number of classes is determined by Sturge's rule). A Mantel test is run over these distance classes versus the ecological distance matrix. The Mantel statistics obtained from each of these tests are then plotted in a correlogram. A filled-in point on the plot indicates that the Mantel statistic was statistically significant (you may provide what alpha to use).

For more information and examples pertaining to this script, please refer to the accompanying tutorial, which can be found at http://qiime.org/tutorials/distance_matrix_comparison.html.



**Usage:** :file:`compare_distance_matrices.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	`-`-method
		Matrix correlation method to use. Valid options: [mantel, partial_mantel, mantel_corr]
	-i, `-`-input_dms
		The input distance matrices, comma-separated. WARNING: Only symmetric, hollow distance matrices may be used as input. Asymmetric distance matrices, such as those obtained by the UniFrac Gain metric (i.e. `beta_diversity.py <./beta_diversity.html>`_ -m unifrac_g), should not be used as input
	-o, `-`-output_dir
		Path to the output directory
	
	**[OPTIONAL]**
		
	-n, `-`-num_permutations
		The number of permutations to perform when calculating the p-value [default: 100]
	-s, `-`-sample_id_map_fp
		Map of original sample ids to new sample ids [default: None]
	-t, `-`-tail_type
		The type of tail test to perform when calculating the p-value. Valid options: [two sided, less, greater] Two sided is a two-tailed test, while less tests for r statistics less than the observed r statistic, and greater tests for r statistics greater than the observed r statistic. Only applies when method is mantel [default: two sided]
	-a, `-`-alpha
		The value of alpha to use when denoting significance in the correlogram plot. Only applies when method is mantel_corr
	-g, `-`-image_type
		The type of image to produce. Valid options: [png, svg, pdf]. Only applies when method is mantel_corr [default: pdf]
	`-`-variable_size_distance_classes
		If this option is supplied, each distance class will have an equal number of distances (i.e. pairwise comparisons), which may result in variable sizes of distance classes (i.e. each distance class may span a different range of distances). If this option is not supplied, each distance class will have the same width, but may contain varying numbers of pairwise distances in each class. This option can help maintain statistical power if there are large differences in the number of distances in each class. See Darcy et al. 2011 (PLoS ONE) for an example of this type of correlogram. Only applies when method is mantel_corr [default: False]
	-c, `-`-control_dm
		The control matrix. Only applies (and is *required*) when method is partial_mantel. [default: None]


**Output:**


Mantel: One file is created containing the Mantel 'r' statistic and p-value.

Partial Mantel: One file is created in the output directory, which contains the partial Mantel statistic and p-value.

Mantel Correlogram: Two files are created in the output directory: a text file containing information about the distance classes, their associated Mantel statistics and p-values, etc. and an image of the correlogram plot.



**Partial Mantel:**

Performs a partial Mantel test on two distance matrices, using a third matrix as a control. Runs 99 permutations to calculate the p-value.

::

	compare_distance_matrices.py --method partial_mantel -i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt -c PH_dm.txt -o partial_mantel_out -n 99

**Mantel:**

Performs a Mantel test on all pairs of four distance matrices, including 999 permutations for each test.

::

	compare_distance_matrices.py --method mantel -i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt,weighted_unifrac_even100_dm.txt,unweighted_unifrac_even100_dm.txt -o mantel_out -n 999

**Mantel Correlogram:**

This example computes a Mantel correlogram on two distance matrices using 999 permutations in each Mantel test. Output is written to the mantel_correlogram_out directory.

::

	compare_distance_matrices.py --method mantel_corr -i unweighted_unifrac_dm.txt,PH_dm.txt -o mantel_correlogram_out -n 999


