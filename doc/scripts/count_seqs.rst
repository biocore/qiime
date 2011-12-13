.. _count_seqs:

.. index:: count_seqs.py

*count_seqs.py* -- 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`count_seqs.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		The input filepaths (comma-separated)
	
	**[OPTIONAL]**
		
	-o, `-`-output_fp
		The output filepath [default: write to stdout]
	`-`-suppress_errors
		Suppress warnings about missing files [default: False]


**Output:**




Count the sequences in a fasta file and write results to stdout.

::

	count_seqs.py -i in.fasta

Count the sequences in a fasta file and a fastq file and write results to file. Note that fastq files can only be processed if they end with .fastq -- all other files are assumed to be fasta.

::

	count_seqs.py -i in1.fasta,in2.fastq -o seq_counts.txt

Count the sequences all .fasta files in current directory and write results to stdout. Note that -i option must be quoted.

::

	count_seqs.py -i "*.fasta"


