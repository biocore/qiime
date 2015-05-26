.. _plot_semivariogram:

.. index:: plot_semivariogram.py

*plot_semivariogram.py* -- Fits a model between two distance matrices and plots the result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Fits a spatial autocorrelation model between two matrices and plots the result. This script will work with two distance matrices but will ignore the 0s at the diagonal and the values that go to N/A. See `distance_matrix_from_mapping.py <./distance_matrix_from_mapping.html>`_.


**Usage:** :file:`plot_semivariogram.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-x, `-`-input_path_x
		Path to distance matrix to be displayed in the x axis
	-y, `-`-input_path_y
		Path to distance matrix to be displayed in the y axis
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation
	
	**[OPTIONAL]**
		
	-b, `-`-binning
		Binning ranges. Format: [increment,top_limit], when top_limit is -1=infinitum; you can specify several ranges using the same format, i.e. [2.5,10][50,-1] will set two bins, one from 0-10 using 2.5 size steps and from 10-inf using 50 size steps. Note that the binning is used to clean the plots (reduce number of points) but ignored to fit the model. [default: None]
	`-`-ignore_missing_samples
		This will overpass the error raised when the matrices have different sizes/samples
	`-`-x_max
		X axis max limit [default: auto]
	`-`-x_min
		X axis min limit [default: auto]
	`-`-y_max
		Y axis max limit [default: auto]
	`-`-y_min
		Y axis min limit [default: auto]
	-X, `-`-x_label
		Label for the x axis [default: Distance Dissimilarity (m)]
	-Y, `-`-y_label
		Label for the y axis [default: Community Dissimilarity]
	-t, `-`-fig_title
		Title of the plot [default: Semivariogram]
	`-`-dot_color
		Dot color for plot, more info: http://matplotlib.sourceforge.net/api/pyplot_api.html [default: white]
	`-`-dot_marker
		Dot color for plot, more info: http://matplotlib.sourceforge.net/api/pyplot_api.html [default: o]
	`-`-line_color
		Line color for plot, more info: http://matplotlib.sourceforge.net/api/pyplot_api.html [default: blue]
	`-`-dot_alpha
		Alpha for dots, more info: http://matplotlib.sourceforge.net/api/pyplot_api.html [default: 1]
	`-`-line_alpha
		Alpha for dots, more info: http://matplotlib.sourceforge.net/api/pyplot_api.html [default: 1]
	`-`-model
		Model to be fitted to the data. Valid choices are:nugget, exponential, gaussian, periodic, linear. [default: exponential]
	-p, `-`-print_model
		Print in the title of the plot the function of the fit. [default: False]
	-c, `-`-category
		Category to color each of the trajectories when you have multiple treatments [default: None]
	-m, `-`-mapping_fp
		Metadata mapping file, only used when coloring by a category, a file with the legends and color coding will be created with the suffix legend [default: None]


**Output:**

The resulting output file consists of a pdf image containing the plot between the two distances matrices and the fitted model


**Fitting:**

For this script, the user supplies two distance matrices (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_), along with the output filename (e.g. semivariogram), and the model to fit, as follows:

::

	plot_semivariogram.py -x distance.txt -y unifrac.txt -o semivariogram_exponential.png

Modify the the default method to gaussian

::

	plot_semivariogram.py -x distance.txt -y unifrac.txt --model gaussian -o semivariogram_gaussian.png

**Color semivariograms by a category in the metadata mapping file:**

Using a header name in the mapping file (Time), create two separate semivariograms in the same plot, an accompanying file with the color coding will be created(categories_legend.eps), both the legends and the plot will be in eps format.

::

	plot_semivariogram.py -y unweighted_unifrac_dm.txt -x time_dm.txt --model gaussian -m Fasting_Map.txt -o categories.eps -c Treatment


