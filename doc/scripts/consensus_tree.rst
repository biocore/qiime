.. _consensus_tree:

.. index:: consensus_tree.py

*consensus_tree.py* -- This script outputs a majority consensus tree given a collection of input trees.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`consensus_tree.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Input folder containing trees
	-o, `-`-output_fname
		The output consensus tree filepath
	
	**[OPTIONAL]**
		
	-s, `-`-strict
		Use only nodes occurring >50% of the time [default: False]


**Output:**

The output is a newick formatted tree compatible with most standard tree viewing programs


**basic usage: given a directory of trees 'jackknifed_trees', compute the majority consensus and save as a newick formatted text file:**

::

	consensus_tree.py -i jackknifed_trees -o consensus_tree.tre


