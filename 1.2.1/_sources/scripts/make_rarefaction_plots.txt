.. _make_rarefaction_plots:

.. index:: make_rarefaction_plots.py

*make_rarefaction_plots.py* -- Generate Rarefaction Plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from `collate_alpha.py <./collate_alpha.html>`_, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied collated alpha-diversity files in a folder or a comma-separated list of files, by passing the "-i" option.  Be aware that this script produces many images for the interactive html pages, so you may choose to not create these pages. The user may also supply optional arguments like an image type (-g), and a resolution (-d).


**Usage:** :file:`make_rarefaction_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Name of folder containing rarefaction files, takes output from `collate_alpha.py <./collate_alpha.html>`_.  The user can also supply a list of files, which are comma-separated. [REQUIRED]
	-m, `-`-map_fname
		Name of mapping file [REQUIRED]
	
	**[OPTIONAL]**
		
	-b, `-`-colorby
		Name of columns to make rarefaction graphs of, comma delimited no spaces.
	-p, `-`-prefs_path
		Preferences file for coloring of columns.
	-k, `-`-background_color
		Background color for graphs.
	-g, `-`-imagetype
		Extension for image type choose from (png, svg, pdf).  WARNING: Some formats may not properly open in your browser! [default: png]
	-d, `-`-resolution
		Output image resolution, [default: 75]
	-y, `-`-ymax
		This is the ymax value to be used for the plots, so you can compare rarefaction plots between two different analyses [default: None]
	-w, `-`-webpage
		This is allows to user to not create the webpage, which may be slow with large datasets [default: True]
	-o, `-`-output_dir
		Path to the output directory


**Output:**

The result of this script produces a folder and within that folder there is a sub-folder containing image files. Within the main folder, there is an html file.


**Default Example:**

For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:

::

	make_rarefaction_plots.py -i collated_alpha/ -m mapping_file.txt

**Specify Image Type and Resolution:**

Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:

::

	make_rarefaction_plots.py -i collated_alpha/ -m mapping_file.txt -d 180 -g pdf

**Use Prefs File:**

You can also supply a preferences file "-p", as follows

::

	make_rarefaction_plots.py -i collated_alpha/ -m mapping_file.txt -d 180 -p prefs.txt

**Set Background Color:**

Alternatively, you can set the plot background "-k", as follows: a preferences file "-p", as follows

::

	make_rarefaction_plots.py -i collated_alpha/ -m mapping_file.txt -k black

**Generate raw data without interactive webpages:**

The user can choose to not create an interactive webpage ("-w" option).  This is for the case, where the user just wants the average plots and the raw average data.

::

	make_rarefaction_plots.py -i collated_alpha/ -m mapping_file.txt -w


