.. _sort_denoiser_output:

.. index:: sort_denoiser_output.py

*sort_denoiser_output.py* -- Sort denoiser output by cluster size.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This scripts is used prior to OTU picking when combining several separately denoised data sets. It sorts the FASTA file by cluster size, such that the OTU pciker now which are the most likely the best OTU centroids.


**Usage:** :file:`sort_denoiser_output.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-input_fasta_fp
		Path to the input fasta file
	-o, `-`-output_file
		The output filename


**Output:**

A standard FASTA file


**Example Usage:**

::

	sort_denoiser_output.py -f denoised_seqs.fasta -o denoised_seqs_sorted.fasta


