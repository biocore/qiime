.. _plot_taxa_summary:

.. index:: plot_taxa_summary.py

*plot_taxa_summary.py* -- Make taxaonomy summary charts based on taxonomy assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script automates the construction of pie, bar and area charts showing the breakdown of taxonomy by given levels. The script creates an html file for each chart type for easy visualization. It uses the taxonomy or category counts from `summarize_taxa.py <./summarize_taxa.html>`_ for combined samples by level (-i) and user specified labels for each file passed in (-l). Output will be written to the user specified folder (-o) the, where the default is the current working directory. The user  can also specify the number of categories displayed for within a single pie chart, where the rest are grouped together as the  'other category' using the (-n) option, default is 20.



**Usage:** :file:`plot_taxa_summary.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-counts_fname
		Input comma-separated list of summarized taxa filepaths (i.e results from `summarize_taxa.py <./summarize_taxa.html>`_) [REQUIRED]
	
	**[OPTIONAL]**
		
	-l, `-`-labels
		Comma-separated list of taxonomic levels (e.g. Phylum,Class,Order)  [default=None]
	-n, `-`-num_categories
		The maximum number of taxonomies to show in each pie chart. All additional taxonomies are grouped into an "other" category. NOTE: this functionality only applies to the pie charts. [default: 20]
	-o, `-`-dir_path
		Output directory
	-b, `-`-colorby
		This is the categories to color by in the plots from the metadata mapping file. The categories must match the name of a  column header in the mapping file exactly and multiple categories can be list by comma separating them without spaces. [default=None]
	-p, `-`-prefs_path
		Input user-generated preferences filepath. NOTE: This is a file with a dictionary containing preferences for the analysis. The key taxonomy_coloring is used for the coloring. [default: None]
	-k, `-`-background_color
		This is the background color to use in the plots (black or white) [default: white]
	-d, `-`-dpi
		This is the resolution of the plot. [default: 80]
	-x, `-`-x_width
		This is the width of the x-axis to use in the plots. [default: 12]
	-y, `-`-y_height
		This is the height of the y-axis to use in the plots. [default: 6]
	-w, `-`-bar_width
		This the width of the bars in the bar graph and should be a number between 0 and 1. NOTE: this only applies to the bar charts. [default: 0.75]
	-t, `-`-type_of_file
		This is the type of image to produce (i.e. pdf,svg,png). [default: pdf]
	-c, `-`-chart_type
		This is the type of chart to plot (i.e. pie, bar or area). The user has the ability to plot multiple types, by using a comma-separated list (e.g. area,pie) [default: area,bar]
	-r, `-`-resize_nth_label
		Make every nth label larger than the other lables. This is for large area and bar charts where the font on the x-axis is small. This requires an integer value greater than 0. [default: 0]
	-s, `-`-include_html_legend
		Include HTML legend. If present, the writing of the legend in the html page is included. [default: False]
	-m, `-`-include_html_counts
		Include HTML counts. If present, the writing of the counts in the html table is included [default: False]
	-a, `-`-label_type
		Label type ("numeric" or "categorical").  If the label type is defined as numeric, the x-axis will be scaled accordingly. Otherwise the x-values will treated categorically and be evenly spaced [default: categorical].


**Output:**

The script generates an output folder, which contains several files. For each pie chart there is a png and a pdf file. The best way to view all of the pie charts is by opening up the file taxonomy_summary_pie_chart.html.


**Examples:**

If you wish to run the code using default parameters, you must supply a counts file (phylum.txt) along with the taxon level label (Phylum), the type(s) of charts to produce, and an output directory, by using the following command:

::

	plot_taxa_summary.py -i phylum.txt -l phylum -c pie,bar,area -o phylum_charts/

If you want to make charts for multiple levels at a time (phylum.txt,class.txt,genus.txt) use the following command:

::

	plot_taxa_summary.py -i phylum.txt,class.txt,genus.txt -l Phylum,Class,Genus -c pie,bar,area -o phylum_class_genus_charts/

Additionally, if you would like to display on a set number of taxa ("-n 10") in the pie charts, you can use the following command:

::

	plot_taxa_summary.py -i class.txt -l Class -c pie -n 10 -o class_pie_n10_charts/

If you would like to display generate pie charts for specific samples, i.e. sample 'PC.636' and sample 'PC.635' that are in the counts file header, you can use the following command:

::

	plot_taxa_summary.py -i class.txt -l Class -b PC.636,PC.635 -o sample_charts/


