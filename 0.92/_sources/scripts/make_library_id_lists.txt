.. _make_library_id_lists:

.. index:: make_library_id_lists.py

*make_library_id_lists.py* -- Make library id lists
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Makes a list of the ids corresponding to each library represented in the input fasta file. Assumes that the libraries are the output of `split_libraries.py <./split_libraries.html>`_ and that they contain the 454 read id for each sequence as is standard in the `split_libraries.py <./split_libraries.html>`_ output. Produces a separate file for each library. These are used to retrieve the corresponding reads from the sff files for SRA deposition.


**Usage:** :file:`make_library_id_lists.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta
		The path to a FASTA file containing input sequences
	
	**[OPTIONAL]**
		
	-s, `-`-screened_rep_seqs
		The path to a FASTA file containing screened representative seqs[DEFAULT: None]
	-u, `-`-otus
		The path to an OTU file mapping OTUs onto rep seqs[DEFAULT: None]
	-o, `-`-outdir
		 The base directory to save results (one file per library).
	-f, `-`-field
		Index of space-delimited field to read id from [DEFAULT: 1]
	`-`-debug
		Show debug output.


**Output:**

This script produces a separate file for each library.


**Example:**

Create a list containing library ids for a fasta file (seqs.fna):

::

	make_library_id_lists.py -i seqs.fna -o results/


