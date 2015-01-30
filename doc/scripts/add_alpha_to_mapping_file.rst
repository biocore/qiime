.. _add_alpha_to_mapping_file:

.. index:: add_alpha_to_mapping_file.py

*add_alpha_to_mapping_file.py* -- Add alpha diversity data to a metadata mapping file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Add alpha diversity data to a mapping file for use with other QIIME scripts, i. e. make_emperor. The resulting mapping file will contain three new columns per metric in the alpha diversity data; the first column being the raw value, the second being a normalized raw value and the third one a label classifying the bin where this value fits based on the normalized value.


**Usage:** :file:`add_alpha_to_mapping_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-alpha_fps
		Alpha diversity data with one or multiple metrics i. e. the output of `alpha_diversity.py <./alpha_diversity.html>`_. This can also be a comma-separated list of collated alpha diversity file paths i. e. the output of `collate_alpha.py <./collate_alpha.html>`_, when using collated alpha diversity data the --depth option is required
	-m, `-`-mapping_fp
		Mapping file to modify by adding the alpha diversity data
	
	**[OPTIONAL]**
		
	-o, `-`-output_mapping_fp
		Filepath for the modified mapping file [default: mapping_file_with_alpha.txt]
	-b, `-`-number_of_bins
		Number of bins [default: 4].
	-x, `-`-missing_value_name
		Bin prefix name for the sample identifiers that exist in the mapping file (mapping_fp) but not in the alpha diversity file (alpha_fp) [default: N/A].
	`-`-binning_method
		Select the method name to create the bins, the options are 'equal' and 'quantile'. Both methods work over the normalized alpha diversity values. On the one hand 'equal' will assign the bins on equally spaced limits, depending on the value of --number_of_bins i. e. if you select 4 the limits will be [0.25, 0.50, 0.75]. On the other hand 'quantile' will select the limits based on the --number_of_bins i. e. the limits will be the quartiles if 4 is selected [default: equal].
	`-`-depth
		Select the rarefaction depth to use when the alpha_fps refers to collated alpha diversity file(s) i. e. the output of `collate_alpha.py <./collate_alpha.html>`_. All the iterations contained at this depth will be averaged to form a single mean value [default: highest depth available].
	`-`-collated_input
		Use to specify that the -i option is composed of collated alpha diversity data.


**Output:**

The result of running this script is a metadata mapping file that will include 3 new columns per alpha diversity metric included in the alpha diversity file. For example, with an alpha diversity file with only PD_whole_tree, the new columns will PD_whole_tree_alpha, PD_whole_tree_normalized and PD_whole_tree_bin.


**Adding alpha diversity data:**

Add the alpha diversity values to a mapping file and classify the normalized values into 4 bins, where the limits will be  0 < x <= 0.25 for the first bin 0.25 < x <= 0.5 for the second bin, 0.5 < x <= 0.75 for the third bin and 0.75 < x <= 1 for the fourth bin.

::

	add_alpha_to_mapping_file.py -i adiv_pd.txt -m mapping.txt -b 4 -o alpha_mapping.txt

**Adding alpha diversity data with the quantile method:**

Add the alpha diversity values to a mapping file and classify the normalized values using the quartiles of the distribution of this values.

::

	add_alpha_to_mapping_file.py -i adiv_pd.txt -m mapping.txt -b 4 -o alpha_mapping_quantile.txt --binning_method=quantile

**Adding collated alpha diversity data:**

Add the mean of the alpha diversity values at a specified rarefaction depth, this case is for use with the output of `collate_alpha.py <./collate_alpha.html>`_. It is recommended that the filenames are the name of the metric used in each file.

::

	add_alpha_to_mapping_file.py -i 'shannon.txt,chao1.txt' -m mapping.txt -b 4 -o collated_alpha_mapping.txt --depth=49 --collated_input


