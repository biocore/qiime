.. _principal_coordinates:

.. index:: principal_coordinates.py

*principal_coordinates.py* -- Principal Coordinates Analysis (PCoA)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Principal Coordinate Analysis (PCoA) is commonly used to compare groups of samples based on phylogenetic or count-based distance metrics (see section on `beta_diversity.py <./beta_diversity.html>`_).


**Usage:** :file:`principal_coordinates.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Path to the input distance matrix file(s) (i.e., the output from `beta_diversity.py <./beta_diversity.html>`_). Is a directory for batch processing and a filename for a single file operation.
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation


**Output:**

The resulting output file consists of the principal coordinate (PC) axes (columns) for each sample (rows). Pairs of PCs can then be graphed to view the relationships between samples. The bottom of the output file contains the eigenvalues and % variation explained for each PC.


**PCoA (Single File):**

For this script, the user supplies a distance matrix (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_), along with the output filename (e.g.  beta_div_coords.txt), as follows:

::

	principal_coordinates.py -i beta_div.txt -o beta_div_coords.txt

**PCoA (Multiple Files):**

The script also functions in batch mode if a folder is supplied as input (e.g. from `beta_diversity.py <./beta_diversity.html>`_ run in batch). No other files should be present in the input folder - only the distance matrix files to be analyzed. This script operates on every distance matrix file in the input directory and creates a corresponding principal coordinates results file in the output directory, e.g.:

::

	principal_coordinates.py -i beta_div_weighted_unifrac/ -o beta_div_weighted_pcoa_results/


