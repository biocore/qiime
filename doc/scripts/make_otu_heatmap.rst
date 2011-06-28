.. _make_otu_heatmap:

.. index:: make_otu_heatmap.py

*make_otu_heatmap.py* -- Make heatmap of OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Once an OTU table has been generated, it can be visualized using a heatmap. In these heatmaps each row corresponds to an OTU, and each column corresponds to a sample. The higher the relative abundance of an OTU in a sample, the more intense the color at the corresponsing position in the heatmap. By default, the OTUs (rows) will be clustered by UPGMA hierarchical clustering, and the samples (columns) will be presented in the order in which they appear in the OTU table. Alternatively, the user may pass in a tree to sort the OTUs (rows) or samples (columns), or both. For samples, the user may also pass in a mapping file. If the user passes in a mapping file and a metadata category, samples (columns in the heatmap) will be grouped by category value and subsequently clustered within each group.


**Usage:** :file:`make_otu_heatmap.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to the output directory
	-t, `-`-otu_tree
		Tree file to be used for sorting OTUs in the heatmap
	-m, `-`-map_fname
		Metadata mapping file to be used for sorting Samples in the heatmap
	-c, `-`-category
		Metadata category for annotating samples, used only in --make_image is specified.
	-s, `-`-sample_tree
		Tree file to be used for sorting samples (e.g, output from `upgma_cluster.py <./upgma_cluster.html>`_). If both this and the sample mapping file are provided, the mapping file is ignored.
	`-`-no_log_transform
		Data will not be log-transformed. Without this option, all zeros will be set to a small value (default is 1/2 the smallest non-zero entry). Data will be translated to be non-negative after log transform, and num_otu_hits will be set to 0.
	`-`-absolute_abundance
		Do not normalize samples to sum to 1.[default False]
	`-`-log_eps
		Small value to replace zeros for log transform. [default: 1/2 the smallest non-zero entry].


**Output:**

The heatmap image is located in the specified output directory. It is formatted as a PDF file.


**Examples:**

Using default values:

::

	make_otu_heatmap.py -i otu_table.txt

Different output directory (i.e., "otu_heatmap"):

::

	make_otu_heatmap.py -i otu_table.txt -o otu_heatmap

Sort the heatmap columns by Sample ID's then you should supply the mapping file, as follows:

::

	make_otu_heatmap.py -i otu_table.txt -o otu_heatmap -m mapping_file.txt

Sort the heatmap columns by Sample ID's and the heatmap rows by the order of tips in the tree, you can supply a tree as follows:

::

	make_otu_heatmap.py -i otu_table.txt -o otu_heatmap -m mapping_file.txt -t tree_file.txt

Group the heatmap columns by metadata category (e.g., GENDER), then cluster within each group:

::

	make_otu_heatmap.py -i otu_table.txt -o otu_heatmap -m mapping_file.txt -c 'GENDER'


