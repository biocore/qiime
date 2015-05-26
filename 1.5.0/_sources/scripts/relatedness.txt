.. _relatedness:

.. index:: relatedness.py

*relatedness.py* -- Calculate NRI and NTI using formulas from Phylocom 4.2/3.41
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script calculates NRI and NTI from a path to a Newick formatted tree and a path to a comma separated list of ids in that tree that form the group whose NRI/NTI you want to test. The tree is not required to have distances. If none are found script will use the number of nodes (self inclusive) as their distance from one another. NRI and NTI are calculated as described in the Phylocom manual, not as in Webb 2002, or Webb 2000. The Phylocom manual is freely available on the web and Webb 2002 can be found in the Annual Review of Ecology and Systematics: Phylogenies and Community Ecology Webb 2002.


**Usage:** :file:`relatedness.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-t, `-`-tree_fp
		The tree filepath
	-g, `-`-taxa_fp
		Taxa list filepath
	
	**[OPTIONAL]**
		
	-i, `-`-iters
		Number of iterations to use for sampling tips without replacement (null model 2 community sampling, see see http://bodegaphylo.wikispot.org/Community_Phylogenetics). [default: 1000]
	-m, `-`-methods
		Comma-separated list of metrics to calculate. [default: nri,nti]


**Output:**

Outputs a value for specified tests


**Calculate both NRI and NTI from the given tree and group of taxa:**

::

	relatedness.py -t gg_tree.tre -i ids.txt -m nri,nti

**Calculate only NRI:**

::

	relatedness.py -t gg_tree.tre -i ids.txt -m nri

**Calculate only NTI using a different number of iterations:**

::

	relatedness.py -t gg_tree.tre -i ids.txt -m nti -i 100


