.. _adjust_seq_orientation:

.. index:: adjust_seq_orientation.py

*adjust_seq_orientation.py* -- Get the reverse complement of all sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Write the reverse complement of all seqs in seqs.fasta (-i) to seqs_rc.fasta (default, change output_fp with -o). Each sequence description line will have ' RC' appended to the end of it (default,
leave sequence description lines untouched by passing -r):


**Usage:** :file:`adjust_seq_orientation.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-o, `-`-output_fp
		The output filepath
	-r, `-`-retain_seq_id
		Leave seq description lines untouched [default: append " RC" to seq description lines]


**Output:**




**Example:**

Reverse complement all sequences in seqs.fna and write result to seqs_rc.fna

::

	adjust_seq_orientation.py -i seqs.fna


