.. _split_fasta_on_sample_ids:

.. index:: split_fasta_on_sample_ids.py

*split_fasta_on_sample_ids.py* -- Split a single post-split_libraries.py fasta file into per-sample fasta files.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Split a single post-`split_libraries.py <./split_libraries.html>`_ fasta file into per-sample fasta files. This script requires that the sequences identitifers are in post-`split_libraries.py <./split_libraries.html>`_ format (i.e., SampleID_SeqID). A fasta file will be created for each unique SampleID.


**Usage:** :file:`split_fasta_on_sample_ids.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		The input fasta file to split
	-o, `-`-output_dir
		The output directory [default: None]
	
	**[OPTIONAL]**
		
	`-`-buffer_size
		The number of sequences to read into memory before writing to file (you usually won't need to change this) [default: 500]


**Output:**

This script will produce an output directory with as many files as samples.


Split seqs.fna into one fasta file per sample and store the resulting fasta files in 'out'

::

	split_fasta_on_sample_ids.py -i seqs.fna -o out/


