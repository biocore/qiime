.. _simsam:

.. index:: simsam.py

*simsam.py* -- Simulate samples for each sample in an OTU table, using a phylogenetic tree.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

 This script makes n samples related to each sample in an input otu table

An input OTU table with 3 samples and n=2 will result in an output otu table with 6 samples total: 3 clusters of 2 related samples.

To simulate each of new samples, this script uses a sample in the input OTU table, and for each OTU in that sample the script
traverses rootward on the tree a distance specified by '-d' to a point x. It then randomly selects a tip that decends from x,
(call that new tip 'o2'), and reassigns all observations of the original OTU to the tip/OTU 'o2'.



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


Make 3 related sample for each sample in otu_table.biom.

::

	simsam.py -i otu_table.biom -t rep_set.tre -o otu_table.simsam.biom -d .001 -n 3


