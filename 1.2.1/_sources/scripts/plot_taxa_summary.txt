.. _plot_taxa_summary:

.. index:: plot_taxa_summary.py

*plot_taxa_summary.py* -- Make taxaonomy summary charts based on taxonomy assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script automates the construction of pie, bar and area charts showing the breakdown of taxonomy by given levels. The script creates an html file for each chart type for easy visualization. It uses the taxonomy or category counts from `summarize_taxa.py <./summarize_taxa.html>`_ for combined samples by level (-i) and user specified labels for each file passed in (-l). Output will be written to the user specified folder (-o) the, where the default is the current working directory. The user can also specify the number of categories displayed for within a single pie chart, where the rest are grouped together as the 'other category' using the (-n) option, default is 20.



**Usage:** :file:`plot_taxa_summary.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_files
		List of files with sample counts by taxonomy [REQUIRED]
	-l, `-`-labels
		List of labels (i.e. Phylum, Class) [REQUIRED]
	
	**[OPTIONAL]**
		
	-n, `-`-num
		Maximum number of individual categories in each pie chart. All additional categories are grouped into an "other" category. NOTE: this is only used for the pie charts. [default: 20]
	-o, `-`-dir-prefix
		Output folder
	-b, `-`-colorby
		This is the samples to make charts for in the counts files from `summarize_taxa.py <./summarize_taxa.html>`_. The sample name must match the name of a sample id in the header of the counts file exactly and multiple categories can be list by comma separating them without spaces.  [default: None]
	-p, `-`-prefs_path
		This is the user-generated preferences file. NOTE: This is a file with a dictionary containing preferences for the analysis. The label taxonomy_coloring is used for the coloring, see example prefs file preferences_file. [default: None]
	-k, `-`-background_color
		This is the background color to use in the plots. [default: None]
	-d, `-`-dpi
		This is the dpi to use in the plots. [default: 80]
	-x, `-`-x_width
		This is the width to use in the plots. [default: 12]
	-y, `-`-y_height
		This is the height to use in the plots. [default: 6]
	-w, `-`-bar_width
		This the width of the bars in the bar graph and should be a number between 0 and 1. NOTE: this is only used for bar charts. [default: 0.75]
	-t, `-`-type_of_file
		This is the filename suffix to use for each high-res plot. (i.e. pdf,svg,png) [default: pdf]
	-c, `-`-chart_type
		Type of chart to plot (i.e. pie, bar or area). The user has the ability to plot multiple types, by using a comma-separated list (e.g. area,pie)  [default: area,bar]
	-r, `-`-resize_nth_label
		This is for large area and bar charts where the font on the x-axis is small. This allows you to set every nth label to be larger on the x-axis.This requires an integer value greater than 0.[default: 0]


**Output:**

The script generates an output folder, which contains several files. For each pie chart there is a png and a pdf file. The best way to view all of the pie charts is by opening up the file taxonomy_summary_pie_chart.html.


**Examples:**

If you wish to run the code using default parameters, you must supply a counts file (Class.txt) along with the taxon level label (Class) and the type(s) of chart, by using the following command:

::

	plot_taxa_summary.py -i Class.txt -l Class -c pie,bar,area

If you want to make charts for multiple levels at a time (phylum.txt,class.txt,genus.txt) use the following command:

::

	plot_taxa_summary.py -i phylum.txt,class.txt,genus.txt -l phylum,class,genus -c pie,bar,area

If you want specify an output directory (e.g. "output_charts/", regardless of whether the directory exists, use the following command:

::

	plot_taxa_summary.py -i Class.txt -l Class -c pie,bar,area -o output_charts/

Additionally, if you would like to display on a set number of taxa ("-n 10") in the pie charts, you can use the following command:

::

	plot_taxa_summary.py -i Class.txt -l Class -c pie -o pie_charts/ -n 10

If you would like to display generate pie charts for samples samples: 'sample1' and 'sample2' that are in the counts file header, you can use the following command:

::

	plot_taxa_summary.py -i Class.txt -l Class -o pie_charts/ -b sample1,sample2


