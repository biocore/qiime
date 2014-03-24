.. _make_fastq:

.. index:: make_fastq.py

*make_fastq.py* -- Make fastq file for ERA submission from paired fasta and qual files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The ERA currently requires a separate fastq file for each library, split by library id. This code takes the output from `split_libraries.py <./split_libraries.html>`_ and the corresponding qual files, pulls the qual info by id, and writes everything either to one file or to per-library files.

The fastq format for each record is as follows:

- @seq_id [and optional description]
- seq as bases + [optionally with repeat of seq_id and repeat line]
- qual scores as string of chr(33+qual)



**Usage:** :file:`make_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-input_fasta_fp
		Path to the input fasta file
	-q, `-`-qual
		Names of qual files, comma-delimited
	
	**[OPTIONAL]**
		
	-o, `-`-result_fp
		Path to store results [default: <input_sequences_filename>.fastq]
	-s, `-`-split
		Make separate file for each library [default:False]


**Output:**

This script creates separate fastq files for each library.


**Example:**

Take input fasta file input_fasta_filepath and qual file input_qual_filepath: make separate file for each library (with the -s option: assumes that the fasta file is the output of `split_libraries.py <./split_libraries.html>`_ or similar script):

::

	make_fasta.py -f input_fasta_filepath -q input_qual_filepath -s


