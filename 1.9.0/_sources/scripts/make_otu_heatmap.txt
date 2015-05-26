.. _make_otu_heatmap:

.. index:: make_otu_heatmap.py

*make_otu_heatmap.py* -- Plot heatmap of OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script visualizes an OTU table as a heatmap where each row corresponds to an OTU and each column corresponds to a sample. The higher the relative abundance of an OTU in a sample, the more intense the color at the corresponsing position in the heatmap. By default, the OTUs (rows) will be clustered by UPGMA hierarchical clustering, and the samples (columns) will be presented in the order in which they appear in the OTU table. Alternatively, the user may supply a tree to sort the OTUs (rows) or samples (columns), or both. The user may also pass in a mapping file for sorting samples. If the user passes in a mapping file and a metadata category, samples (columns) will be grouped by category value and subsequently clustered within each group.


**Usage:** :file:`make_otu_heatmap.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-t, `-`-otu_tree
		Tree file to be used for sorting OTUs in the heatmap
	-m, `-`-map_fname
		Metadata mapping file to be used for sorting Samples in the heatmap.
	-c, `-`-category
		Metadata category for sorting samples. Samples will be clustered within each category level using euclidean UPGMA.
	-s, `-`-sample_tree
		Tree file to be used for sorting samples (e.g, output from `upgma_cluster.py <./upgma_cluster.html>`_). If both this and the sample mapping file are provided, the mapping file is ignored.
	-g, `-`-imagetype
		Type of image to produce (i.e. png, svg, pdf) [default: pdf]
	`-`-no_log_transform
		Data will not be log-transformed. Without this option, all zeros will be set to a small value (default is 1/2 the smallest non-zero entry). Data will be translated to be non-negative after log transform, and num_otu_hits will be set to 0.
	`-`-suppress_row_clustering
		No UPGMA clustering of OTUs (rows) is performed. If --otu_tree is provided, this flag is ignored.
	`-`-suppress_column_clustering
		No UPGMA clustering of Samples (columns) is performed. If --map_fname is provided, this flag is ignored.
	`-`-absolute_abundance
		Do not normalize samples to sum to 1 [default: False]
	`-`-color_scheme
		Color scheme for figure. see http://matplotlib.org/examples/color/colormaps_reference.html for choices [default: YlGn]
	`-`-width
		Width of the figure in inches [default: 5]
	`-`-height
		Height of the figure in inches [default: 5]
	`-`-dpi
		Resolution of the figure in dots per inch [default: value of savefig.dpi in matplotlibrc file]
	`-`-obs_md_category
		Observation metadata category to plot [default: taxonomy]
	`-`-obs_md_level
		The level of observation metadata to plot for hierarchical metadata [default: lowest level]


**Output:**

A single output file is created containing the heatmap of the OTU table (a PDF file by default).


Generate a heatmap as a PDF using all default values:

::

	make_otu_heatmap.py -i otu_table.biom -o heatmap.pdf

Generate a heatmap as a PNG:

::

	make_otu_heatmap.py -i otu_table.biom -o heatmap.png -g png

Sort the heatmap columns (samples) by the order of samples in the mapping file

::

	make_otu_heatmap.py -i otu_table.biom -o heatmap_sorted_samples.pdf -m mapping_file.txt

Sort the heatmap columns (samples) by the order of samples in the mapping file, and sort the heatmap rows by the order of tips in the tree:

::

	make_otu_heatmap.py -i otu_table.biom -o heatmap_sorted.pdf -m mapping_file.txt -t rep_set.tre

Group the heatmap columns (samples) by metadata category (e.g., Treatment), then cluster within each group:

::

	make_otu_heatmap.py -i otu_table.biom -o heatmap_grouped_by_Treatment.pdf -m mapping_file.txt -c Treatment


