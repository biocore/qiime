.. _detrend:

.. index:: detrend.py

*detrend.py* -- Detrend Principal Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Ordination plots (e.g. principal coordinates analysis) of samples that lay along a naturally occurring gradient (e.g. depth, time, pH) often exhibit a curved shape known as the "arch" or "horseshoe" effect. This can cause samples near the endpoints of the gradient to appear closer to one another than would be expected. This script will attempt to remove any (compounded) quadratic curvature in a set of 2D coordinates. If requested, it will also report an evaluation of the association of the transformed coordinates with a known gradient.


**Usage:** :file:`detrend.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		Path to read PCoA/PCA/ordination table
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to output directory [default: .]
	-m, `-`-map_fp
		Path to metadata file [default: None]
	-c, `-`-gradient_variable
		Column header for gradient variable in metadata table [default: None]
	-r, `-`-suppress_prerotate
		Suppress pre-rotation of the coordinates for optimal detrending; not pre-rotating assumes that the curvature is symmetrical across the vertical axis [default: False]


**Output:**

Output is detrended PCoA matrices.


**Examples:**

The simplest usage takes as input only a table of principal coordinates:

::

	detrend.py -i $PWD/pcoa.txt -o detrending

One may also include a metadata file with a known real-valued gradient as one of the columns. In this case, the output folder will include a text file providing a summary of how well the analysis fit with the hypothesis that the primary variation is due to the gradient (in this case, "DEPTH"):

::

	detrend.py -i $PWD/pcoa.txt -m map.txt -c DEPTH -o detrending

Note that if you provide a real-valued known gradient the script will prerotate the first two axes of the PCoA coords in order to achieve optimal alignment with that gradient. This can be disabled with "-r":

::

	detrend.py -i $PWD/pcoa.txt -m map.txt -c DEPTH -o detrending -r


