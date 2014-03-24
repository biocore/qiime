.. _nmds:

.. index:: nmds.py

*nmds.py* -- Nonmetric Multidimensional Scaling (NMDS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Nonmetric Multidimensional Scaling (NMDS) is commonly used to compare groups of samples based on phylogenetic or count-based distance metrics (see section on `beta_diversity.py <./beta_diversity.html>`_).


**Usage:** :file:`nmds.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Path to the input distance matrix file(s) (i.e., the output from `beta_diversity.py <./beta_diversity.html>`_). Is a directory for batch processing and a filename for a single file operation.
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation
	
	**[OPTIONAL]**
		
	-d, `-`-dimensions
		Number of dimensions of NMDS spacedefault: 2


**Output:**

The resulting output file consists of the NMDS axes (columns) for each sample (rows). Pairs of NMDS axes can then be graphed to view the relationships between samples. The bottom of the output file contains the stress of the ordination.


**NMDS (Single File):**

For this script, the user supplies a distance matrix (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_), along with the output filename (e.g. beta_div_coords.txt), as follows:

::

	nmds.py -i beta_div.txt -o beta_div_coords.txt

**NMDS (Multiple Files):**

The script also functions in batch mode if a folder is supplied as input (e.g. from `beta_diversity.py <./beta_diversity.html>`_ run in batch). No other files should be present in the input folder - only the distance matrix files to be analyzed. This script operates on every distance matrix file in the input directory and creates a corresponding nmds results file in the output directory, e.g.:

::

	nmds.py -i beta_div_weighted_unifrac/ -o beta_div_weighted_nmds_results/


