.. _make_prefs_file:

.. index:: make_prefs_file.py

*make_prefs_file.py* -- Generate preferences file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates a preferences (prefs) file, which can be passed to `make_2d_plots.py <./make_2d_plots.html>`_ and `make_3d_plots.py <./make_3d_plots.html>`_. The prefs file allows for defining the monte_carlo distance, gradient coloring of continuous values in the 2D and 3D plots, the ball size scale for all the samples and the color of the arrow and the line of the arrow for the procrustes analysis. Currently there is only one color gradient: red to blue.


**Usage:** :file:`make_prefs_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map_fname
		This is the metadata mapping file [default=None]
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-b, `-`-mapping_headers_to_use
		Mapping fields to use in prefs file [default: ALL]
	-k, `-`-background_color
		This is the backgroundcolor to  use in the plots. [default: black]
	-d, `-`-monte_carlo_dists
		Monte carlo distanceto use for each sample header [default: 10]
	-i, `-`-input_taxa_file
		Summarized taxa file with samplecounts by taxonomy (resulting file from `summarize_taxa.py <./summarize_taxa.html>`_)
	-s, `-`-ball_scale
		Scale factor for the size of each ball in the plots [default: 1.0]
	-l, `-`-arrow_line_color
		Arrow line color forprocrustes analysis. [default: white]
	-a, `-`-arrow_head_color
		Arrow head color forprocrustes analysis. [default: red]


**Output:**

The result of this script is a text file, containing coloring preferences to be used by `make_2d_plots.py <./make_2d_plots.html>`_ and `make_3d_plots.py <./make_3d_plots.html>`_.


**Examples:**

To make a prefs file, the user is required to pass in a user-generated mapping file using "-m" and an output filepath, using "-o". When using the defaults, the script will use ALL categories from the mapping file, set the background to black and the monte_carlo distances to 10.

::

	make_prefs_file.py -m mapping.txt -o prefs_out.txt

If the user would like to use specified categories ('SampleID,Individual') or combinations of categories ('SampleID&&Individual'), they will need to use the -b option, where each category is comma delimited, as follows:

::

	make_prefs_file.py -m mapping.txt -b "SampleID,Treatment,SampleID&&Treatment" -o prefs_out_1.txt

If the user would like to change the background color for their plots, they can pass the '-k' option, where the colors: black and white can be used for 3D plots and many additional colors can be used for the 2D plots, such as cyan, pink, yellow, etc.: 

::

	make_prefs_file.py -m mapping.txt -k white -o prefs_out_white.txt

If the user would like to change the monte_carlo distances, they can pass the '-d' option as follows: 

::

	make_prefs_file.py -m mapping.txt -d 15 -o prefs_out_d15.txt

If the user would like to add a list of taxons they can pass the '-i' option, which is the resulting taxa file from `summarize_taxa.py <./summarize_taxa.html>`_, as follows: 

::

	make_prefs_file.py -m mapping.txt -i taxa_level_3.txt -o prefs_out_taxa_l3.txt

If the user would like to add the ball size scale they can pass the '-s' option as follows: 

::

	make_prefs_file.py -m mapping.txt -s 3 -o prefs_out_s3.txt

If the user would like to add the head and line color for the arrows in the procrustes analysis plot they can pass the '-a' and '-l' options as follows: 

::

	make_prefs_file.py -m mapping.txt -a black -l blue -o prefs_out_procrustes.txt


