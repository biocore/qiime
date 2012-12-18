.. _compare_3d_plots:

.. index:: compare_3d_plots.py

*compare_3d_plots.py* -- Plot several PCoA files on the same 3D plot
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates a 3D plot comparing two or more sets of principal coordinates using as input two or more principal coordinates files. Edges are drawn in the plot connecting samples with the same ID across different principal coordinates files. The user can also include a file listing the edges to be drawn in the plot, in which case the user may submit any number of principal coordinates files (including one). If the user includes the edges file, the sample IDs need not match between principal coordinates files.

The principal_coordinates coordinates files are obtained by applying "`principal_coordinates.py <./principal_coordinates.html>`_" to a file containing beta diversity measures. The beta diversity files are optained by applying "`beta_diversity.py <./beta_diversity.html>`_" to an OTU table. One may apply "`transform_coordinate_matrices.py <./transform_coordinate_matrices.html>`_" to the principal_coordinates coordinates files before using this script to compare them.


**Usage:** :file:`compare_3d_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-coord_fnames
		This is comma-separated list of the paths to the principal coordinates files (i.e., resulting file from `principal_coordinates.py <./principal_coordinates.html>`_), e.g 'pcoa1.txt,pcoa2.txt'
	-m, `-`-map_fname
		This is the user-generated mapping file [default=None]
	
	**[OPTIONAL]**
		
	-b, `-`-colorby
		This is a list of the categories to color by in the plots from the user-generated mapping file. The categories must match the name of a column header in the mapping file exactly and multiple categories can be list by comma separating them without spaces. The user can also combine columns in the mapping file by separating the categories by "&&" without spaces [default=None]
	-a, `-`-custom_axes
		This is a category or list of categories from the user-generated mapping file to use as a custom axis in the plot.  For instance, if there is a pH category and one would like to see the samples plotted on that axis instead of PC1, PC2, etc., one can use this option.  It is also useful for plotting time-series data [default: None]
	-p, `-`-prefs_path
		This is the user-generated preferences file. NOTE: This is a file with a dictionary containing preferences for the analysis. See `make_prefs_file.py <./make_prefs_file.html>`_. [default: None]
	-k, `-`-background_color
		This is the background color to use in the plots (Options are 'black' or 'white'. [default: None]
	-e, `-`-edges_file
		A file where each line contains two sample IDs separated by a whitespace character; for each pair of sample IDs, an edge will be drawn from the first sample to the second sample. [default: None]
	`-`-serial
		Connect the 1st set of points to the 2nd, the 2nd to the 3rd, etc. Default behavior is to connect each set of points back to the 1st set. This flag is ignored if the user specifies an edges file.
	-o, `-`-output_dir
		Path to the output directory


**Output:**

This script results in a folder containing an html file which displays the 3D Plots generated.


**Example 1:**

Compare two pca/pcoa files in the same 3d plot where each sample ID is assigned its own color:

::

	compare_3d_plots.py -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -m $PWD/Fasting_Map.txt

**Example 2:**

Compare two pca/pcoa files in the same 3d plot with two coloring schemes (Treatment and DOB):

::

	compare_3d_plots.py -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -m $PWD/Fasting_Map.txt -b 'Treatment,DOB'

**Example 3:**

Compare two pca/pcoa files in the same 3d plot for a combination of label headers from a mapping file: 

::

	compare_3d_plots.py -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -m $PWD/Fasting_Map.txt -b 'Treatment&&DOB' -o $PWD/test/

**Example 4:**

Pass in a list of desired edges and only one pca/pcoa file: 

::

	compare_3d_plots.py -i $PWD/raw_pca_data1.txt -e $PWD/edges.txt -m Fasting_Map.txt -b 'Treatment&&DOB' -o $PWD/test2/

**Example 5:**

Pass in a list of desired edges and only one pca/pcoa file: 

::

	compare_3d_plots.py -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -e $PWD/edges.txt -m $PWD/Fasting_Map.txt -b 'Treatment&&DOB' -o $PWD/test3/


