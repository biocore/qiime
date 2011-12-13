.. _make_distance_boxplots:

.. index:: make_distance_boxplots.py

*make_distance_boxplots.py* -- Creates boxplots to compare distances                                     between categories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script creates boxplots that allow for the comparison between different categories found within the mapping file. The boxplots that are created show the distances within all samples of a field value, as well as between different field values. Individual within and individual between distances are also plotted.

For more information and examples pertaining to this script, please refer to the accompanying tutorial, which can be found at http://qiime.org/tutorials/creating_distance_comparison_plots.html.



**Usage:** :file:`make_distance_boxplots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		The mapping filepath
	-o, `-`-output_dir
		Path to the output directory
	-d, `-`-distance_matrix_fp
		Input distance matrix filepath (i.e. the result of `beta_diversity.py <./beta_diversity.html>`_)
	-f, `-`-fields
		Comma-separated list of fields to compare, where the list of fields should be in quotes (e.g. "Field1,Field2,Field3")
	
	**[OPTIONAL]**
		
	-g, `-`-imagetype
		Type of image to produce (i.e. png, svg, pdf) [default: pdf]
	`-`-save_raw_data
		Store raw data used to create boxplots in tab-delimited files [default: False]
	`-`-suppress_all_within
		Suppress plotting of "all within" boxplot [default: False]
	`-`-suppress_all_between
		Suppress plotting of "all between" boxplot [default: False]
	`-`-suppress_individual_within
		Suppress plotting of individual "within" boxplot(s) [default: False]
	`-`-suppress_individual_between
		Suppress plotting of individual "between" boxplot(s) [default: False]
	`-`-y_min
		The minimum y-axis value in the resulting plot. If "auto", it is automatically calculated [default: 0]
	`-`-y_max
		The maximum y-axis value in the resulting plot. If "auto", it is automatically calculated [default: 1]
	`-`-width
		Width of the output image in inches. If not provided, a "best guess" width will be used [default: auto]
	`-`-height
		Height of the output image in inches [default: 6]
	`-`-transparent
		Make output images transparent (useful for overlaying an image on top of a colored background) [default: False]
	`-`-whisker_length
		Length of the whiskers as a function of the IQR. For example, if 1.5, the whiskers extend to 1.5 * IQR. Anything outside of that range is seen as an outlier [default: 1.5]
	`-`-box_width
		Width of each box in plot units [default: 0.5]
	`-`-box_color
		The color of the boxes. Can be any valid matplotlib color string, such as "black", "magenta", "blue", etc. See http://matplotlib.sourceforge.net/api/colors_api.html for more examples of valid color strings that may be used [default: same as plot background, which is white unless --transparent is enabled]
	`-`-sort
		Sort boxplots by increasing median. If no sorting is applied, boxplots will be grouped logically as follows: all within, all between, individual within, and individual between [default: False]


**Output:**

Images of the plots are written to                                     the specified output directory (one image                                     per field). The raw data used in the                                     plots can optionally be written into                                     tab-delimited files (one file per field).


**Compare distances between Fast and Control for Treatment field:**

This example will generate an image with boxplots for all within and all between distances for the field Treatment, and will also include plots for individual within (e.g. Control vs. Control, Fast vs. Fast) and individual between (e.g. Control vs. Fast). The generated plot image PDF will be written to the output directory 'out_files'.

::

	make_distance_boxplots.py -d dist_matrix.txt -m map.txt -f "Treatment" -o out_files

**Only plot individual distances:**

This example will generate a PNG of all of individual distances (within and between) for the Treatment field.

::

	make_distance_boxplots.py -d dist_matrix.txt -m map.txt -f "Treatment" -o out_files -g png --suppress_all_within --suppress_all_between

**Save raw data:**

This example will generate an SVG image of the boxplots and also output the plotting data to a tab-delimited file.

::

	make_distance_boxplots.py -d dist_matrix.txt -m map.txt -f "Treatment" -o out_files -g svg --save_raw_data


