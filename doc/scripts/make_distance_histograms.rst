.. _make_distance_histograms:

.. index:: make_distance_histograms.py

*make_distance_histograms.py* -- Make distance histograms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


To visualize the distance between samples and/or categories in the metadata
mapping file, the user can generate histograms to represent the distances
between samples. This script generates an HTML file, where the user can compare
the distances between samples based on the different categories associated to
each sample in the metadata mapping file.

Distance histograms provide a way to compare different categories and see which
tend to have larger/smaller distances than others. For example, in a hand
study, you may want to compare the distances between hands to the distances
between individuals (with the file "hand_distances.txt" using the parameter -d
hand_distances.txt). The categories are defined in the metadata mapping file
(specified using the parameter -m hand_map.txt). If you want to look at the
distances between hands and individuals, choose the "Hand" field and
"Individual" field (using the parameter --fields Hand,Individual (notice the
fields are comma-delimited)). For each of these groups of distances a
histogram is made. The output is an HTML file which is created in the
"Distance_Histograms" directory (using the parameter -o Distance_Histograms to
specify output directory) where you can look at all the distance histograms
individually, and compare them between each other.



**Usage:** :file:`make_distance_histograms.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-d, `-`-distance_matrix_file
		Input distance matrix filepath (i.e. the result of `beta_diversity.py <./beta_diversity.html>`_). WARNING: Only symmetric, hollow distance matrices may be used as input. Asymmetric distance matrices, such as those obtained by the UniFrac Gain metric (i.e. `beta_diversity.py <./beta_diversity.html>`_ -m unifrac_g), should not be used as input
	-m, `-`-map_fname
		Input metadata mapping filepath.
	
	**[OPTIONAL]**
		
	-p, `-`-prefs_path
		Input user-generated preferences filepath. NOTE: This is a file with a dictionary containing preferences for the analysis. This dictionary must have a "Fields" key mapping to a list of desired fields. [default: None]
	-o, `-`-dir_path
		Output directory. [default: ./]
	-k, `-`-background_color
		Background color for use in the plots (black or white) [default: white]
	`-`-monte_carlo
		Deprecated: pass --monte_carlo_iters > 0 to enable
	`-`-suppress_html_output
		Suppress HTML output. [default: False]
	-f, `-`-fields
		Comma-separated list of fields to compare, where the list of fields should be in quotes (e.g. "Field1,Field2,Field3"). Note: if this option is passed on the command-line, it will overwrite the fields in prefs file. [default: first field in mapping file is used]
	`-`-monte_carlo_iters
		Number of iterations to perform for Monte Carlo analysis. [default: 0; No monte carlo simulation performed]


**Output:**


The result of this script will be a folder containing images and/or an HTML
file (with appropriate javascript files), depending on the user-defined
parameters.



**Distance histograms example:**

In the following command, the user supplies a distance matrix (i.e. the resulting file from `beta_diversity.py <./beta_diversity.html>`_), the user-generated metadata mapping file and one category "Treatment" to generate distance histograms.

::

	make_distance_histograms.py -d unweighted_unifrac_dm.txt -m Fasting_Map.txt --fields Treatment -o example1

**Multiple categories:**

For comparison of multiple categories (e.g. Treatment, DOB), you can use the following command (separating each category with a comma).

::

	make_distance_histograms.py -d unweighted_unifrac_dm.txt -m Fasting_Map.txt --fields Treatment,DOB -o example2

**Suppress HTML output:**

By default, HTML output is automatically generated. If the user would like to suppress the HTML output, you can use the following command.

::

	make_distance_histograms.py -d unweighted_unifrac_dm.txt -m Fasting_Map.txt --fields Treatment --suppress_html_output -o example3

**Preferences file:**

You can provide your own preferences file (prefs.txt) with the following command. If a preferences file is supplied, you do not need to supply fields on the command-line.

::

	make_distance_histograms.py -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -p prefs.txt -o example4


