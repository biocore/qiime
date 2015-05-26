.. _make_distance_comparison_plots:

.. index:: make_distance_comparison_plots.py

*make_distance_comparison_plots.py* -- Creates plots comparing distances between sample groupings
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script creates plots (bar charts, scatter plots, or box plots) that
allow for the comparison between samples grouped at different field states
of a mapping file field.

This script can work with any field in the mapping file, and it can compare
any number of field states to all other field states within that field.
This script may be especially useful for fields that represent a time series,
because a plot can be generated showing the distances between samples at
certain timepoints against all other timepoints.

For example, a time field might contain the values 1, 2, 3, 4, and 5, which
label samples that are from day 1, day 2, day 3, and so on. This time field
can be specified when the script is run, as well as the timepoint(s) to
compare to every other timepoint. For example, two comparison groups
might be timepoints 1 and 2. The resulting plot would contain timepoints for
days 3, 4, and 5 along the x-axis, and at each of those timepoints, the
distances between day 1 and that timepoint would be plotted, as well as the
distances between day 2 and the timepoint.

The script also performs two-sample t-tests for all pairs of distributions to
help determine which distributions are significantly different from each other.

Tip: the script tries its best to fit everything into the plot, but there are
cases where plot elements may get cut off (e.g. if axis labels are extremely
long), or things may appear squashed, cluttered, or too small (e.g. if
there are many boxplots in one plot). Increasing the width and/or height of the
plot (using --width and --height) usually fixes these problems.

For more information and examples pertaining to this script, please refer to
the accompanying tutorial, which can be found at
http://qiime.org/tutorials/creating_distance_comparison_plots.html.



**Usage:** :file:`make_distance_comparison_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		The mapping filepath
	-o, `-`-output_dir
		Path to the output directory
	-d, `-`-distance_matrix_fp
		Input distance matrix filepath (i.e. the result of `beta_diversity.py <./beta_diversity.html>`_). WARNING: Only symmetric, hollow distance matrices may be used as input. Asymmetric distance matrices, such as those obtained by the UniFrac Gain metric (i.e. `beta_diversity.py <./beta_diversity.html>`_ -m unifrac_g), should not be used as input
	-f, `-`-field
		Field in the mapping file to make comparisons on
	-c, `-`-comparison_groups
		Comma-separated list of field states to compare to every other field state, where the list of field states should be in quotes (e.g. "FieldState1,FieldState2,FieldState3")
	
	**[OPTIONAL]**
		
	-t, `-`-plot_type
		Type of plot to produce ("bar" is bar chart, "scatter" is scatter plot, and "box" is box plot) [default: bar]
	-g, `-`-imagetype
		Type of image to produce (i.e. png, svg, pdf) [default: pdf]
	`-`-save_raw_data
		Store raw data used to create plot in a tab-delimited file [default: False]
	`-`-suppress_significance_tests
		Suppress performing signifance tests between each pair of distributions [default: False]
	-n, `-`-num_permutations
		The number of Monte Carlo permutations to perform when calculating the nonparametric p-value in the significance tests. Must be an integer greater than or equal to zero. If zero, the nonparametric p-value will not be calculated and will instead be reported as "N/A". This option has no effect if --suppress_significance_tests is supplied [default: 0]
	`-`-tail_type
		The type of tail test to compute when calculating the p-values in the significance tests. "high" specifies a one-tailed test for values greater than the observed t statistic, while "low" specifies a one-tailed test for values less than the observed t statistic. "two-sided" specifies a two-tailed test for values greater in magnitude than the observed t statistic. This option has no effect if --suppress_significance_tests is supplied. Valid choices: low or high or two-sided [default: two-sided]
	`-`-width
		Width of the output image in inches [default: 12]
	`-`-height
		Height of the output image in inches [default: 6]
	`-`-x_tick_labels_orientation
		Type of orientation for x-axis tick labels [default: vertical]
	-a, `-`-label_type
		Label type ("numeric" or "categorical"). If the label type is defined as numeric, the x-axis will be scaled accordingly. Otherwise the x-values will treated categorically and will be evenly spaced [default: categorical].
	`-`-y_min
		The minimum y-axis value in the resulting plot. If "auto", it is automatically calculated [default: 0]
	`-`-y_max
		The maximum y-axis value in the resulting plot. If "auto", it is automatically calculated [default: 1]
	`-`-transparent
		Make output images transparent (useful for overlaying an image on top of a colored background ) [default: False]
	`-`-whisker_length
		If --plot_type is "box", determines the length of the whiskers as a function of the IQR. For example, if 1.5, the whiskers extend to 1.5 * IQR. Anything outside of that range is seen as an outlier. If --plot_type is not "box", this option is ignored [default: 1.5]
	`-`-error_bar_type
		If --plot_type is "bar", determines the type of error bars to use. "stdv" is standard deviation and "sem" is the standard error of the mean. If --plot_type is not "bar", this option is ignored [default: stdv]
	`-`-distribution_width
		Width (in plot units) of each individual distribution (e.g. each bar if the plot type is a bar chart, or the width of each box if the plot type is a boxplot) [default: auto]


**Output:**


An image of the plot is written to the specified output directory. The raw data
used in the plots and the results of significance tests can optionally be
written into tab-delimited files that are most easily viewed in a spreadsheet
program such as Microsoft Excel.



**Compare distances between Native and Input samples for each timepoint in the Time field:**

This example will generate a PDF containing a bar chart with the distances between Native samples and every other timepoint, as well as the distances between Input samples and every other timepoint. The output image will be put in the 'out1' directory. For more details about this example input data, please refer to the accompanying tutorial.

::

	make_distance_comparison_plots.py -d forearm_only_unweighted_unifrac_dm.txt -m costello_timeseries_map.txt -f TIME_SINCE_TRANSPLANT -c "Native,Input" -o out1


