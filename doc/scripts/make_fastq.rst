.. _make_fastq:

.. index:: make_fastq.py

*make_fastq.py* -- Make FASTQ file for ERA submission from paired FASTA and QUAL files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The ERA currently requires a separate FASTQ file for each library, split by library id. This code takes the output from `split_libraries.py <./split_libraries.html>`_ and the corresponding QUAL files and produces ERA-compatible FASTQ files.


**Usage:** :file:`make_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-input_fasta_fp
		Path to the input fasta file
	-q, `-`-qual
		Names of QUAL files, comma-delimited
	
	**[OPTIONAL]**
		
	-o, `-`-result_fp
		Path to store results [default: <input_sequences_filename>.fastq]
	-s, `-`-split
		Make separate file for each library [default:False]


**Output:**

Matches QUAL info to FASTA entries by id, and writes FASTQ output to one file or to per-library files.

The FASTQ format for each record is as follows:

@seq_id [and optional description]
seq as bases
+ [and optionally with repeat of seq_id and repeat line]
qual scores as string of chr(33+qual)



**Example:**

Take input FASTA file input_fasta_filepath and QUAL file input_qual_filepath: make separate file for each library (with the -s option: assumes that the FASTA file is the output of `split_libraries.py <./split_libraries.html>`_ or similar script):

::

	make_fastq.py -f $PWD/seqs.fna -q $PWD/Fasting_Example.qual -s


