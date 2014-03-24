.. _make_3d_plot_prefs_file:

.. index:: make_3d_plot_prefs_file.py

*make_3d_plot_prefs_file.py* -- Generate preferences file for 3D plots using Metadata Fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates a preferences (prefs) file, which can be passed via -p to `make_3d_plots.py <./make_3d_plots.html>`_. The prefs file allows for gradient coloring of continuous values in the 3D plots. Currently there is only one color gradient: red to blue.


**Usage:** :file:`make_3d_plot_prefs_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-b, `-`-color_by
		Mapping fields to color by [default: None]
	-p, `-`-output_prefs_fp
		Path to store output file [default: None]


**Output:**

The result of this script is a text file, containing coloring preferences to be used by `make_3d_plots.py <./make_3d_plots.html>`_.


**Example:**

To make a prefs file to be used by `make_3d_plots.py <./make_3d_plots.html>`_, the -b value should be passed in as the same as that passed in via -b to `make_3d_plots.py <./make_3d_plots.html>`_. For multiple fields, the command should be a delimited list of fields. For example the -b string could be "#SampleID,Individual" and output to the file "prefs_out.txt" using the -p parameter.

::

	make_3d_plot_prefs_file.py -b "#SampleID,Individual" -p prefs_out.txt


