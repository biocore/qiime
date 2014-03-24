.. _make_rarefaction_averages:

.. index:: make_rarefaction_averages.py

*make_rarefaction_averages.py* -- Generate Rarefaction Averages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Once the batch alpha diversity files have been collated, you may want to compare the diversity. Using the results from `collate_alpha.py <./collate_alpha.html>`_, you can average rarefaction values across sample metadata in the mapping file using this script.

This script creates a directory of average rarefaction series based on the supplied mapping file (-m) and the supplied rarefaction files (-r) from `collate_alpha.py <./collate_alpha.html>`_.


**Usage:** :file:`make_rarefaction_averages.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map
		Name of mapping file [REQUIRED]
	-r, `-`-rarefaction
		Name of rarefaction file, takes output from collate_alpha OR tab delimited data from a previous run of this script. If using raw data from a previous run, set -x flag. [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-prefs
		Name of columns to make rarefaction graphs of, comma delimited no spaces. Use 'ALL' command to make graphs of all metadata columns. [default: ALL]
	-o, `-`-dir_path
		Directory prefix for all analyses [default: .]


**Output:**

The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be text files of averages for each of the categories specified by the user.


**Default Example:**

For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:

::

	make_rarefaction_averages.py -m Mapping_file.txt -r chao1.txt

**Multiple File Example:**

If you would like to generate plots for multiple files, you can use the following command:

::

	make_rarefaction_averages.py -m Mapping_file.txt -r chao1.txt,PD_whole_tree.txt

**Category Specific Example:**

In the case that you want to make plots for a specific category (i.e., pH), you can use the following command:

::

	make_rarefaction_averages.py -m Mapping_file.txt -r chao1.txt -p pH


