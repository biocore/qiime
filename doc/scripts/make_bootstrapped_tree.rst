.. _make_bootstrapped_tree:

.. index:: make_bootstrapped_tree.py

*make_bootstrapped_tree.py* -- Make bootstrapped tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a tree and bootstrap support file and creates a pdf, colored by bootstrap support.


**Usage:** :file:`make_bootstrapped_tree.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-master_tree
		This is the path to the master tree
	-s, `-`-support
		This is the path to the bootstrap support file
	-o, `-`-output_file
		This is the filename where the output should be written.


**Output:**

The result of this script is a pdf file.


**Example:**

In this example, the user supplies a tree file and a text file containing the jackknife support information, which results in a pdf file:

::

	make_bootstrapped_tree.py -m master_tree.tre -s jackknife_support.txt -o jackknife_samples.pdf


