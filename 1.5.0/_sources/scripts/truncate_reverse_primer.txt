.. _truncate_reverse_primer:

.. index:: truncate_reverse_primer.py

*truncate_reverse_primer.py* -- Takes a demultiplexed fasta file, finds a specified reverse primer sequence, and truncates this primer and subsequent sequences following the reverse primer.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Takes fasta sequences have already been demultiplexed (via `split_libraries.py <./split_libraries.html>`_, `denoise_wrapper.py <./denoise_wrapper.html>`_, `ampliconnoise.py <./ampliconnoise.html>`_, etc.), and have fasta labels that are QIIME format, i.e., SampleID_#, this script will use the SampleID and a mapping file with a ReversePrimer column to find the reverse primer by local alignment and remove this and any subsequent sequence in a filtered output fasta file.


**Usage:** :file:`truncate_reverse_primer.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fasta_fp
		Fasta file.  Needs to have fasta labels in proper demultiplexed format 
	-m, `-`-mapping_fp
		Mapping filepath.  ReversePrimer field required.  Reverse primers need to be in 5'->3' orientation.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Output directory.  Will be created if does not exist.  [default: .]
	-z, `-`-truncate_option
		Truncation option.  The default option, "truncate_only" will try to find the reverse primer to truncate, and if not found, will write the sequence unchanged.  If set to "truncate_remove", sequences where the reverse primer is not found will not be written. [default: truncate_only]
	-M, `-`-primer_mismatches
		Number of mismatches allowed in the reverse primer. [default: 2]


**Output:**

Truncated version of the input fasta file (based on input name with '_reverse_primer_removed' appended) will be generated in the output directory, along with a .log file.


**Example:**

Find, truncate reverse primers from the fasta file seqs.fna, write output fasta file to the reverse_primer_removed directory:

::

	truncate_reverse_primer.py -f seqs.fna -m Fasting_Map.txt -o reverse_primer_removed/


