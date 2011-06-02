.. _make_per_library_sff:

.. index:: make_per_library_sff.py

*make_per_library_sff.py* -- Make per-library sff files from ID lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates per-library sff files using a directory of text files, one per library, which list read ID's to be included.

The ID list files should contain one read ID per line. If a line contains multiple words (separated by whitespace), then only the first word is used. A '>' character is stripped from the beginning of the line, if present. Blank lines in the file are skipped.



**Usage:** :file:`make_per_library_sff.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_sff
		Input sff file (separate multiple files w/ comma)
	-l, `-`-libdir
		Directory containing ID list text files, one per library
	
	**[OPTIONAL]**
		
	-p, `-`-sfffile_path
		Path to sfffile binary [default: use sfffile in $PATH]
	`-`-use_sfftools
		Use external sfffile program instead of equivalent Python routines.
	`-`-debug
		Print debugging output to stdout [default: False]


**Output:**

The result of this script generates sff files for each library.


**Example:**

Make per-library sff files using input.sff and a directory of libs where each file in the directory contains the id lists for each library:

::

	make_per_library_sff.py -i input.sff -l libs


