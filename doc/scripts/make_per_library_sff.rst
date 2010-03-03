.. _make_per_library_sff:

.. index:: make_per_library_sff.py

*make_per_library_sff.py* -- Make per-library sff files from id lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates per-library sff files using the id lists.


**Usage:** :file:`make_per_library_sff.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_sff
		The path to an input sff file (or files: separate w/ comma, no spaces)
	-l, `-`-libdir
		 The directory containing per-library id files
	
	**[OPTIONAL]**
		
	-p, `-`-sfffile_path
		 Path to sfffile binary [default: sfffile]
	`-`-debug
		Print command-line for debugging [default: False]


**Output:**

The result of this script generates sff files for each library.


**Example:**

Make per-library sff files using input.sff and a directory of libs where each file in the directory contains the id lists for each library:

::

	make_per_library_sff.py -i input.sff -l libs


