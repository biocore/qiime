.. _categorized_dist_scatterplot:

.. index:: categorized_dist_scatterplot.py

*categorized_dist_scatterplot.py* -- makes a figure representing average distances between samples, broken down by categories. I call it a 'categorized distance scatterplot'
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

makes a figure representing average distances between samples, broken down by categories. I call it a 'categorized distance scatterplot'. See script usage for more details. The mapping file specifies the relavent data - if you have e.g. 'N/A' values or samples you don't want included, first use `filter_by_metadata.py <./filter_by_metadata.html>`_ to remove unwanted samples from the mapping file, and thus the analysis. Note that the resulting plot will include only samples in both the mapping file AND the distance matrix.


**Usage:** :file:`categorized_dist_scatterplot.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map
		Mapping file
	-d, `-`-distance_matrix
		Distance matrix
	-p, `-`-primary_state
		Samples matching this state will be plotted. E.g.: AgeCategory:Child . See qiime's `filter_by_metadata.py <./filter_by_metadata.html>`_ for more syntax options
	-a, `-`-axis_category
		This will form the horizontal axis of the figure, e.g.: AgeYears . Must be numbers
	-o, `-`-output_path
		Output figure, filename extention determines format. E.g.: "fig1.png" or similar. A "fig1.txt" or similar will also be created with the data underlying the figure
	
	**[OPTIONAL]**
		
	-c, `-`-colorby
		Samples will first be separated by this column of the mapping file. They will be colored by this column of the mapping file, and all comparisons will be done only among samples with the same value in this column. e.g.: Country. You may omit -c, and the samples will not be separated
	-s, `-`-secondary_state
		All samples matching the primary state will be compared to samples matcthing this secondary state. E.g.: AgeCategory:Adult


**Output:**

a figure and the text dat for that figure 


**Canonical Example:**

Split samples by country. Within each country compare each child to all adults. Plot the average distance from that child to all adults, vs. the age of that child

::

	python categorized_dist_scatterplot.py -m map.txt -d unifrac_distance.txt -c Country -p AgeCategory:Child -s AgeCategory:Adult -a AgeYears -o fig1.png

**Example 2:**

Same as above, but compares Child with all other categories (e.g.: NA, Infant, etc.)

::

	python categorized_dist_scatterplot.py -m map.txt -d unifrac_distance.txt -c Country -p AgeCategory:Child -a AgeYears -o fig1.svg


