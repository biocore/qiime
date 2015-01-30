.. _split_sequence_file_on_sample_ids:

.. index:: split_sequence_file_on_sample_ids.py

*split_sequence_file_on_sample_ids.py* -- Split a single post-split_libraries.py fasta (or post-split_libraries_fastq.py fastq) file into per-sample files.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Split a single post-`split_libraries.py <./split_libraries.html>`_ fasta (or post-`split_libraries_fastq.py <./split_libraries_fastq.html>`_ fastq) file into per-sample fasta files. This script requires that the sequences identitifers are in post-`split_libraries.py <./split_libraries.html>`_ format (i.e., SampleID_SeqID). A file will be created for each unique SampleID.


**Usage:** :file:`split_sequence_file_on_sample_ids.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_seqs_fp
		The input fasta file to split
	-o, `-`-output_dir
		The output directory [default: None]
	
	**[OPTIONAL]**
		
	`-`-buffer_size
		The number of sequences to read into memory before writing to file (you usually won't need to change this) [default: 500]
	`-`-file_type
		Type of file. Either fasta or fastq


**Output:**

This script will produce an output directory with as many files as samples.


Split seqs.fna into one fasta file per sample and store the resulting fasta files in 'out'

::

	split_sequence_file_on_sample_ids.py -i seqs.fna -o out/

Split seqs.fastq into one fastq file per sample and store the resulting fastq files in 'out_fastq'

::

	split_sequence_file_on_sample_ids.py -i seqs.fastq --file_type fastq -o out_fastq/


