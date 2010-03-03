.. _make_rarefaction_plots:

.. index:: make_rarefaction_plots.py

*make_rarefaction_plots.py* -- Generate Rarefaction Plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from `collate_alpha.py <./collate_alpha.html>`_, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied mapping file (-m) and the supplied rarefaction files (-r) from `collate_alpha.py <./collate_alpha.html>`_. The user may also supply optional arguments that will only create plots for supplied metadata columns from the mapping file (-p), an image type (-i), and a resolution (-d). If the user would like to suppress html output they can pass the -n flag, and output raw data with the -w flag. The -y option allows the user to supply a maximum value for the yaxis of the plots.


**Usage:** :file:`make_rarefaction_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map
		Name of mapping file [REQUIRED]
	-r, `-`-rarefaction
		Name of rarefaction file, takesoutput from collate_alpha OR tab delimited data from a previous run of this script. If using raw data from a previous run, set -x flag. [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-prefs
		Name of columns to make rarefaction graphs of, comma delimited no spaces. Use 'ALL' command to make graphs of all metadata columns. [default: ALL]
	-n, `-`-no_html
		Suppress html output. [default: False]
	-i, `-`-imagetype
		Extension for image type choose from (jpg, gif, png, svg, pdf). [default: png]
	-d, `-`-resolution
		Output image resolution, [default: 75]
	-o, `-`-dir_path
		Directory prefix for all analyses [default: .]
	-y, `-`-ymax
		Maximum value for y axis, [default: 0] the default value will tell the script to calculate a y axis maximum depending on the data
	-x, `-`-raw_data
		Read in tab delimited, rarefaction graphing data to be plotted.
	-w, `-`-write_raw_data
		Print out tab delimited, rarefaction graphing data that can be read back in and plotted. [default: False]


**Output:**

The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be images for each of the categories specified by the user.


**Default Example:**

For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:

::

	make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt

**Multiple File Example:**

If you would like to generate plots for multiple files, you can use the following command:

::

	make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt,PD_whole_tree.txt

**Category Specific Example:**

In the case that you want to make plots for a specific category (i.e., pH), you can use the following command:

::

	make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt -p

**Specify Image Type and Resolution:**

Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:

::

	make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt -p pH -d 180 -i pdf


