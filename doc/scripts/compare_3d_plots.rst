.. _compare_3d_plots:

.. index:: compare_3d_plots.py

*compare_3d_plots.py* -- Plot two PCoA files on the same 3D plot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates a 3D plot comparing two sets of principle coordinates using as input two principle coordinates files. The principle coordinates files are obtained by applying "`principle_coordinates.py <./principle_coordinates.html>`_" to a file containing beta diversity measures. The beta diversity files are optained by applying "`beta_diversity.py <./beta_diversity.html>`_" to an OTU table. One may apply "`transform_coordinate_matrices.py <./transform_coordinate_matrices.html>`_" to the principle coordinates files before using this script to compare them.


**Usage:** :file:`compare_3d_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-coord_fnames
		This is comma-separated list of the paths to the principal coordinates files (i.e., resulting file from `principal_coordinates.py <./principal_coordinates.html>`_), e.g 'pcoa1.txt,pcoa2.txt'
	
	**[OPTIONAL]**
		
	-m, `-`-map_fname
		This is the user-generated mapping file [default=None]
	-b, `-`-colorby
		This is a list of the categories to color by in the plots from the user-generated mapping file. The categories must match the name of a column header in the mapping file exactly and multiple categories can be list by comma separating them without spaces. The user can also combine columns in the mapping file by separating the categories by "&&" without spaces [default=None]
	-a, `-`-custom_axes
		This is a category or list of categories from the user-generated mapping file to use as a custom axis in the plot.  For instance, if there is a pH category and one would like to see the samples plotted on that axis instead of PC1, PC2, etc., one can use this option.  It is also useful for plotting time-series data [default: None]
	-p, `-`-prefs_path
		This is the user-generated preferences file. NOTE: This is a file with a dictionary containing preferences for the analysis. See `make_3d_plot_prefs_file.py <./make_3d_plot_prefs_file.html>`_. [default: None]
	-o, `-`-dir_path
		This is the location where the resulting output should be written [default=]


**Output:**

This script results in a folder containing an html file which displays the 3D Plots generated.


**Example 1:**

Compare two pca/pcoa files in the same 3d plot where each sample ID is assigned its own color:

::

	compare_3d_plots.py -i 'raw_pca_data1.txt,raw_pca_data2.txt'

**Example 2:**

Compare two pca/pcoa files in the same 3d plot with two coloring schemes (Day and Type):

::

	compare_3d_plots.py -i 'raw_pca_data1.txt,raw_pca_data2.txt' -m input_map.txt -b 'Day,Type'

**Example 3:**

Compare two pca/pcoa files in the same 3d plot for a combination of label headers from a mapping file: 

::

	compare_3d_plots.py -i 'raw_pca_data1.txt,raw_pca_data2.txt' -m input_map.txt -b 'Type&&Day' -o ./test/


