.. _compare_trajectories:

.. index:: compare_trajectories.py

*compare_trajectories.py* -- Run analysis of volatility using a variety of algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script mainly allows performing analysis of volatility on time series data, but they can be applied to any data that contains a gradient. The methods available are RMS (either using 'avg' or 'trajectory'); or the first difference (using 'diff'), or 'wdiff' for a modified first difference algorithm. The trajectories are computed as follows. For 'avg' it calculates the average point within a group and then computes the norm of the distance of each sample from the averaged value. For 'trajectory' each component of the result trajectory is computed as taking the sorted list of samples in the group and taking the norm of the coordinates of the 2nd samples minus the 1st sample, 3rd sample minus 2nd sample and so on. For 'diff' it calculates the norm for all the time-points and then calculates the first difference for each resulting point. For 'wdiff', it calculates the norm for all the time-points and substracts the mean of the next number of elements, specified using the '--window_size' parameters, and the current element.


**Usage:** :file:`compare_trajectories.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		Input ordination results filepath
	-m, `-`-map_fp
		Input metadata mapping filepath
	-c, `-`-categories
		Comma-separated list of category names of the mapping file to use to create the trajectories
	-o, `-`-output_dir
		Name of the output directory to save the results
	
	**[OPTIONAL]**
		
	-s, `-`-sort_by
		Category name of the mapping file to use to sort
	`-`-algorithm
		The algorithm to use. Available methods: ['avg', 'trajectory', 'diff', 'wdiff']. [Default: avg]
	`-`-axes
		The number of axes to account while doing the trajectory specific calculations. We suggest using 3 because those are the ones being displayed in the plots but you could use any number between 1 and number of samples - 1. To use all of them pass 0. [default: 3]
	-w, `-`-weight_by_vector
		Use -w when you want the output to be weighted by the space between samples in the --sort_by column, i. e. days between samples [default: False]
	`-`-window_size
		Use --window_size, when selecting the modified first difference ('wdiff') option for --algorithm. This integer determines the number of elements to be averaged per element subtraction, the resulting trajectory. [default: None]


**Output:**

This script generates two files in the output directory, 'trajectories.txt' and 'trajectories_raw_values.txt'. The 'trajectories.txt' file includes the resulting statistics and a list of categories that did not passed the tests to run the analysis. The 'trajectories_raw_values.txt' file includes the resulting trajectory for each used category.


**Average method:**

Execute the analysis of volatility using the average method, grouping the samples using the 'Treatment' category

::

	compare_trajectories.py -i pcoa_res.txt -m map.txt -c 'Treatment' -o avg_output

**Trajectory method:**

Execute the analysis of volatility using the trajectory method, grouping the samples using the 'Treatment' category and sorting them using the 'time' category

::

	compare_trajectories.py -i pcoa_res.txt -m map.txt -c 'Treatment' --algorithm trajectory -o trajectory_output -s time

**First difference method:**

Execute the analysis of volatility using the first difference method, grouping the samples using the 'Treatment' category, sorting them using the 'time' category and calculating the trajectory using the first four axes

::

	compare_trajectories.py -i pcoa_res.txt -m map.txt -c 'Treatment' --algorithm diff -o diff_output -s time --axes 4

**Window difference method:**

Execute the analysis of volatility using the window difference method, grouping the samples using the 'Treatment' category, sorting them using the 'time' category, weighting the output by the space between samples in the 'time' category and using a window size of three.

::

	compare_trajectories.py -i pcoa_res.txt -m map.txt -c 'Treatment' --algorithm wdiff -o wdiff_output -s time --window_size 3 -w


