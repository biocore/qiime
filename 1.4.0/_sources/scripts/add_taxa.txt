.. _add_taxa:

.. index:: add_taxa.py

*add_taxa.py* -- Add taxa to OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script adds taxa to an OTU table.


**Usage:** :file:`add_taxa.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_file
		Path to read otu file
	-t, `-`-taxonomy_file
		Path to read taxonomy file
	-o, `-`-output_file
		Path to write
	
	**[OPTIONAL]**
		
	-m, `-`-id_map_file
		Path to read seq id to otu map file [default: None]


**Output:**

The result of this script is written to the specified file.


**Example:**

Add taxa to otu file from otus.txt from file taxa.txt:

::

	add_taxa.py -i otus.txt -t taxa.txt


