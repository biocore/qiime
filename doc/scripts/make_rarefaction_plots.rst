.. _make_rarefaction_plots:

.. index:: make_rarefaction_plots.py

*make_rarefaction_plots.py* -- Generate Rarefaction Plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from `make_rarefaction_averages.py <./make_rarefaction_averages.html>`_, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied rarefaction files in the folder given (-i) from `make_rarefaction_averages.py <./make_rarefaction_averages.html>`_. The user may also supply optional arguments like an image type (-i), and a resolution (-d).


**Usage:** :file:`make_rarefaction_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		name of folder containing rarefaction files, takes output from `make_rarefaction_data.py <./make_rarefaction_data.html>`_ [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-rarefactionAve
		name of overall average rarefaction file, takes output from `make_rarefaction_data.py <./make_rarefaction_data.html>`_
	-p, `-`-imagetype
		extension for image type choose from (jpg, gif, png, svg, pdf). [default: png]
	-d, `-`-resolution
		output image resolution, [default: 75]
	-o, `-`-dir_path
		directory prefix for all analyses [default: .]
	-y, `-`-ymax
		maximum value for y axis, [default: 0] the default value will tell the script to calculate a y axis maximum depending on the data


**Output:**

The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be images for each of the categories specified by the user.


**Default Example:**

For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:

::

	make_rarefaction_plots.py -r chao1/

**Specify Image Type and Resolution:**

Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:

::

	make_rarefaction_plots.py -i chao1/ -d 180 -p pdf


