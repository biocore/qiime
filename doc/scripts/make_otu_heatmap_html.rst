.. _make_otu_heatmap_html:

.. index:: make_otu_heatmap_html.py

*make_otu_heatmap_html.py* -- Make heatmap of OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Create an interactive OTU heatmap from an OTU table. This script parses the OTU count table and filters the table by counts per otu (user-specified), then converts the table into a javascript array, which can be loaded into a web application. The OTU heatmap displays raw OTU counts per sample, where the counts are colored based on the contribution of each OTU to the total OTU count present in that sample (blue: contributes low percentage of OTUs to sample; red: contributes high percentage of OTUs). This web application allows the user to filter the otu table by number of counts per otu. The user also has the ability to view the table based on taxonomy assignment. Additional features include: the ability to drag rows (up and down) by clicking and dragging on the row headers; and the ability to zoom in on parts of the heatmap by clicking on the counts within the heatmap.


**Usage:** :file:`make_otu_heatmap_html.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	-o, `-`-output_dir
		Path to the output directory
	
	**[OPTIONAL]**
		
	-n, `-`-num_otu_hits
		Only include OTUs with at least this many sequences. [default: 5]
	-t, `-`-tree
		Path to newick tree where OTUs are tips, used for sorting OTUs in the heatmap
	-m, `-`-map_fname
		Input metadata mapping filepath, used for sorting samples in the heatmap
	`-`-sample_tree
		Path to newick tree where samples are tips (e.g, output from `upgma_cluster.py <./upgma_cluster.html>`_) used for sorting samples in the heatmap. If both this and the metadata mapping file are provided, the mapping file will be ignored.
	`-`-log_transform
		Log-transform the data. All zeros will be set to a small value (default is 1/2 of the smallest non-zero entry). Data will be translated to be non-negative after log transform and the num_otu_hits will be set to 0.
	`-`-log_eps
		Small value to replace zeros when performing log transformation. [default: 1/2 the smallest non-zero entry].


**Output:**

The interactive heatmap is located in OUTPUT_DIR/otu_table.html where OUTPUT_DIR is specified as -o. Safari is recommended for viewing the OTU Heatmap, since the HTML table generation is much faster than Firefox (as of this writing).


**Generate an OTU heatmap:**

By using the default values ("-n 5), you can then use the code as follows:

::

	make_otu_heatmap_html.py -i otu_table.biom -o heatmap/

**Generate a filtered OTU heatmap:**

If you would like to filter the OTU table by a different number of counts per OTU (i.e., 10):

::

	make_otu_heatmap_html.py -i otu_table.biom -n 10 -o heatmap_mc10/

**Generate a sample-sorted OTU heatmap:**

If you would like to sort the heatmap by Sample IDs then you should supply the mapping file, as follows:

::

	make_otu_heatmap_html.py -i otu_table.biom -o heatmap_sample_sorted -m Fasting_Map.txt

**Generate a sample and OTU-sorted OTU heatmap:**

If you would like to sort the heatmap by Sample IDs and the tips in the tree, you can supply a tree as follows:

::

	make_otu_heatmap_html.py -i otu_table.biom -o heatmap_sample_otu_sorted -m Fasting_Map.txt -t rep_set.tre


