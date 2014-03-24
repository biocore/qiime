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
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	-o, `-`-output_fp
		The output filepath


**Output:**

The resulting output file consists of each component (columns) along with the loading for each sample (rows). Pairs of components can then be graphed to view the relationships between samples. The bottom of the output file contains the eigenvalues and % variation explained for each component.


**Example:**

For this script, the user supplies a distance matrix (i.e., resulting file from `beta_diversity.py <./beta_diversity.html>`_), along with the output filename (e.g. beta_div_coords.txt), as follows:

::

	principal_coordinates.py -i beta_div.txt -o beta_div_coords.txt


