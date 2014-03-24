.. _clean_raxml_parsimony_tree:

.. index:: clean_raxml_parsimony_tree.py

*clean_raxml_parsimony_tree.py* -- Remove duplicate tips from Raxml Tree
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script allows the user to remove specific duplicate tips from a Raxml tree.


**Usage:** :file:`clean_raxml_parsimony_tree.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_tree
		The input raxml parsimony tree
	-t, `-`-tips_to_keep
		The input tips to score and retain (comma-separated list)
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-s, `-`-scoring_method
		The scoring method either depth or numtips [default: depth]


**Output:**




**Example (depth):**

For this case the user can pass in input Raxml tree, duplicate tips, and define an output filepath. When using the depth option, only the deepest replicate is kept. 

::

	 clean_raxml_parsimony_tree.py -i raxml_v730_final_placement.tre -t 6 -o raxml_v730_final_placement_depth.tre

**Example (numtips):**

For this case the user can pass in input Raxml tree, duplicate tips, and define an output filepath. When using the numtips option, the replicate with the fewest siblings is kept. 

::

	 clean_raxml_parsimony_tree.py -i raxml_v730_final_placement.tre -t 6 -o raxml_v730_final_placement_numtips.tre -s numtips


