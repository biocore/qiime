.. _simsam:

.. index:: simsam.py

*simsam.py* -- Simulate samples for each sample in an OTU table, using a phylogenetic tree.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`simsam.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table
		The input otu table
	-t, `-`-tree_file
		Tree file
	-o, `-`-output_file
		The output file
	-d, `-`-dissim
		Dissimilarity between nodes up the tree
	-n, `-`-num
		Number of simulated samples per input sample


**Output:**

an otu table, samples are named: 'original_sample_0, original_sample_1 ...'


Make 3 related sample for each sample in in_otu_table.txt.

::

	simsam.py -i in_otu_table.txt -t rep_set.tre -o out_otu_table.txt -d .001 -n 3


