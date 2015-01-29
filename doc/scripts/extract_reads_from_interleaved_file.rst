.. _extract_reads_from_interleaved_file:

.. index:: extract_reads_from_interleaved_file.py

*extract_reads_from_interleaved_file.py* -- Extract reads from an interleaved file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes an interleaved file, like the ones produced by JGI, and outputs a forward and reverse fastq file with the corresponding reads in each file. 


**Usage:** :file:`extract_reads_from_interleaved_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		Path to input forward reads in FASTQ format.
	-o, `-`-output_dir
		Directory to store result files
	
	**[OPTIONAL]**
		
	`-`-forward_read_identifier
		This is the string identifying the forward reads. [default: 1:N:0].
	`-`-reverse_read_identifier
		This is the string identifying the reverse reads. [default: 2:N:0].


**Output:**

A new folder with two fastq files: forward_reads.fastq and reverse_reads.fastq


**Extract reads from an interleaved file:**

::

	 extract_reads_from_interleaved_file.py -i $PWD/reads_to_extract.fastq -o $PWD/extracted_reads


