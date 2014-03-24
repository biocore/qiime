.. _make_otu_heatmap_html:

.. index:: make_otu_heatmap_html.py

*make_otu_heatmap_html.py* -- Make heatmap of OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Once the OTU table has been generated, the user can create an interactive OTU heatmap. This script parses the OTU count table and filters the table by counts per otu (user-specified), then converts the table into a javascript array, which can be loaded into a web application. The OTU heatmap displays raw OTU counts per sample, where the counts are colored based on the contribution of each OTU to the total OTU count present in that sample (blue: contributes low percentage of OTUs to sample; red: contributes high percentage of OTUs). This web application allows the user to filter the otu table by number of counts per otu. The user also has the ability to view the table based on taxonomy assignment. Additional features include: the ability to drag rows (up and down) by clicking and dragging on the row headers; and the ability to zoom in on parts of the heatmap by clicking on the counts within the heatmap.


**Usage:** :file:`make_otu_heatmap_html.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to the output directory
	-n, `-`-num_otu_hits
		This is the minimum number of Samples that an OTU is present in, for an OTU to be kept in the OTU table [default: 5]


**Output:**

The interactive heatmap is located in a randomly generated folder where the name of the folder starts with "otu_table". The resulting folder contains the interactive heatmap (html file) along with a javascript library folder. This web application has been tested in Mozilla Firefox and Safari. Safari is recommended for viewing the OTU Heatmap, since the HTML table generation is much faster.


**Examples:**

By using the default values ("-n 5), you can then use the code as follows:

::

	make_otu_heatmap_html.py -i otu_table.txt

If you would like to filter the OTU table by a different number of counts per OTU (i.e., 10), you can use the following code:

::

	make_otu_heatmap_html.py -i otu_table.txt -n 10

If you would like to specify a different output directory (i.e., "otu_heatmap"), you can use the following code:

::

	make_otu_heatmap_html.py -i otu_table.txt -o otu_heatmap


