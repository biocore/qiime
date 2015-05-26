.. _plot_semivariogram:

.. index:: plot_semivariogram.py

*plot_semivariogram.py* -- Fits a model between two distance matrices and plots the result
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Fits a spatial autocorrelation model between two matrices and plots the result. This script will work with two distance matrices but will ignore the 0s at the diagonal and the values that go to N/A


**Usage:** :file:`plot_semivariogram.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-x, `-`-input_path_x
		Path to distance matrix to be displayed in the x axis
	-y, `-`-input_path_y
		Path to distance matrix to be displayed in the y axis
	-m, `-`-model
		Model to be fitted to the data. Valid choices are:nugget, exponential, gaussian, periodic. [default: exponential]
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation
	
	**[OPTIONAL]**
		
	-b, `-`-binning
		Binning ranges. Format: [increment,top_limit], when top_limit is -1=infinitum; you can specify several ranges using the same format, i.e. [2.5,10][50,-1] will set two bins, one from 0-10 using 2.5 size steps and from 10-inf using 50 size steps. Note that the binning is used to clean the plots (reduce number of points) but ignored to fit the model. [default: None]


**Output:**

The resulting output file consists of a pdf image containing the plot between the two distances matrices and the fitted model


**Fitting:**

For this script, the user supplies two distance matrices (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_), along with the output filename (e.g. semivariogram), and the model to fit, as follows:

::

	plot_semivariogram.py -x distance.txt -y unifrac.txt -m exponential -o semivariogram.png


