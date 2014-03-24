.. _plot_rank_abundance_graph:

.. index:: plot_rank_abundance_graph.py

*plot_rank_abundance_graph.py* -- plot rank-abundance curve
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Plot a set of rank-abundance graphs from an OTU table and a set of sample names. Multiple graphs will be plotted into the same figure, in order to allow for an easy comparison across samples.


**Usage:** :file:`plot_rank_abundance_graph.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	-s, `-`-sample_name
		Name of the sample to plot. Use "*" to plot all.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Name of output directory. [default: random]
	-a, `-`-absolute_counts
		Plot absolute abundance values instead of relative [default: False]
	-n, `-`-no-legend
		Do not draw a legend [default: False]
	-x, `-`-x_linear_scale
		Draw x axis in linear scale [default: False]
	-y, `-`-y_linear_scale
		Draw y axis in linear scale [default: False]
	-f, `-`-file_type
		Save plot using this image type. Choice of pdf, svg, png, eps [default: pdf]


**Output:**




**Single graph example:**

Plot the rank-abundance curve of one sample using a linear scale for the x_axis:

::

	plot_rank_abundance_graph.py -i otu_table.txt  -s 'Sample1' -x -v

**multiple graph example:**

Plot the rank-abundance curve of several sample:

::

	plot_rank_abundance_graph.py -i otu_table.txt  -s 'Sample1,Sample3,Sample5' -x  -v

**multiple graph example:**

Plot the rank-abundance curve of all samples in an OTU table:

::

	plot_rank_abundance_graph.py -i otu_table.txt  -s '*' -x -f eps -v


