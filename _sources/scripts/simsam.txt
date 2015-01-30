.. _simsam:

.. index:: simsam.py

*simsam.py* -- Simulate samples for each sample in an OTU table, using a phylogenetic tree.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script makes n samples related to each sample in an input otu table 
An input OTU table with 3 samples and n=2 will result in an output OTU table with 6 samples total: 3 clusters of 2 related samples.
To simulate each of the new samples, this script uses a sample in the input OTU table, and for each OTU in that sample the script traverses rootward on the tree a distance specified by '-d' to a point x. It then randomly selects a tip that decends from x, (call that new tip 'o2'), and reassigns all observations of the original OTU to the tip/OTU 'o2'.


**Usage:** :file:`simsam.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table
		The input otu table
	-t, `-`-tree_file
		Tree file
	-o, `-`-output_dir
		Path to the output directory
	-d, `-`-dissim
		Dissimilarity between nodes up the tree, as a single value or comma-separated list of values
	-n, `-`-num
		Number of simulated samples per input sample, as a single value or comma-separated list of values
	
	**[OPTIONAL]**
		
	-m, `-`-mapping_fp
		The mapping filepath. If provided, an output mapping file containing the replicated sample IDs (with all other metadata columns copied over) will also be created [default: None]


**Output:**


The output directory will contain an OTU table with samples named:
'original_sample_0, original_sample_1 ...'

If a mapping file is provided via -m, an output mapping file containing the
replicated sample IDs (with all other metadata columns copied over) will also
be created.



Create an OTU table with 3 related samples for each sample in otu_table.biom with dissimilarities of 0.001.

::

	simsam.py -i otu_table.biom -t rep_set.tre -o simsam_out1 -d .001 -n 3

Create OTU tables with 2, 3 and 4 related samples for each sample in otu_table.biom with dissimilarities of 0.001 and 0.01. Additionally create new mapping files with metadata for each of the new samples for use in downstream analyses.

::

	simsam.py -i otu_table.biom -t rep_set.tre -o simsam_out2 -d .001,.01 -n 2,3,4 -m map.txt


