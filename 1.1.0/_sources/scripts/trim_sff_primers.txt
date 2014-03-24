.. _trim_sff_primers:

.. index:: trim_sff_primers.py

*trim_sff_primers.py* -- Trim sff primers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Finds the technical read regions for each library, and resets the left trim.


**Usage:** :file:`trim_sff_primers.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-l, `-`-libdir
		The directory containing per-library sff files
	-m, `-`-input_map
		Path to the input mapping file describing the libraries
	
	**[OPTIONAL]**
		
	-p, `-`-sfffile_path
		Path to sfffile binary [default: sfffile]
	-q, `-`-sffinfo_path
		Path to sffinfo binary [default: sffinfo]
	`-`-debug
		Print command-line output for debugging [default: False]


**Output:**

This script replaces the original sff files with the trimmed versions.


**Simple example:**

Trim a directory of per-sff files in sff_dir (-l sff_dir/) using an input map (-m input_map.txt). This script uses the sff utility binaries which must be in your path.

::

	trim_sff_primers.py -l sff_dir/ -m input_map.txt


