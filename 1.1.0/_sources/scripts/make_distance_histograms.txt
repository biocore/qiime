.. _make_distance_histograms:

.. index:: make_distance_histograms.py

*make_distance_histograms.py* -- Make distance histograms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To visualize the distance between samples and/or categories in the metadata mapping file, the user can generate histograms to represent the distances between samples. This script generates an HTML file, where the user can compare the distances between samples based on the different categories associated to each sample in the metadata mapping file. 


**Usage:** :file:`make_distance_histograms.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-d, `-`-distance_matrix_file
		Path to distance matrix file.
	-m, `-`-map_fname
		This is the metadata mapping file  [default=None]
	
	**[OPTIONAL]**
		
	-p, `-`-prefs_path
		This is the user-generated preferences file. NOTE: This is a file with a dictionary containing preferences for the analysis.  This dict must have a "Fields" key mapping to a list of desired fields. [default: None]
	-o, `-`-dir_path
		Directory to output data for all analyses. [default: .]
	-k, `-`-background_color
		This is the     background color to use in the plots (Options are 'black' or 'white'.     [default: white]
	`-`-monte_carlo
		Perform Monte Carlo analysis on distances.  [Default: False]
	`-`-html_output
		Write output in HTML format. [Default: False]
	-f, `-`-fields
		Comma delimited list of fields to compare.  This overwrites fields in prefs file.  If this is not provided, the first field in metadata mapping file will be used.  Usage: --fields Field1,Field2,Field3
	`-`-monte_carlo_iters
		Number of iterations to perform for Monte Carlo analysis. [default: 10]


**Output:**

The result of this script will be a folder containing images and/or an html file (with appropriate javascript files), depending on the user-defined parameters.


**Examples:**

Distance Histograms are a way to compare different categories and see which tend to have larger/smaller distances than others. For example, in the hand study, you may want to compare the distances between hands to the distances between individuals (with the file "hand_distances.txt" using the parameter -d hand_distances.txt). The categories are defined in the metadata mapping file (specified using the parameter -m hand_map.txt). If you want to look at the distances between hands and individuals, choose the "Hand" field and "Individual" field (using the parameter --fields Hand,Individual (notice the fields are comma delimited)). For each of these groups of distances a histogram is made. The output is a HTML file ("QIIME_Distance_Histograms.html" when the parameter --html_output is specified) which is created in the "Distance_Histograms" directory (using the parameter -o Distance_Histograms to specify output directory) where you can look at all the distance histograms individually, and compare them between each other.

In the following command, the user only supplies a distance matrix (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_), the user-generated metadata mapping file and one category (e.g. pH):

::

	make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH

For comparison of multiple categories (e.g. pH, salinity), you can use the following command:

::

	make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH,salinity

If the user would like to write the result to a dynamic HTML, you can use the following command:

::

	make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH --html_output

In the case that the user generates their own preferences file (prefs.txt), they can use the following command:

::

	make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt -p prefs.txt

Note: In the case that a preferences file is passed, the user does not need to supply fields in the command-line.


