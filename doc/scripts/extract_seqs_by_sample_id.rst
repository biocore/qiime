.. _extract_seqs_by_sample_id:

.. index:: extract_seqs_by_sample_id.py

*extract_seqs_by_sample_id.py* -- Extract sequences based on the SampleID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script creates a fasta file which will contain only sequences that ARE associated with the supplied sampleIDs(-s), OR all sequences that are NOT associated with the supplied sampleIDs (-n)


**Usage:** :file:`extract_seqs_by_sample_id.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-s, `-`-sample_ids
		Comma-separated sample_ids to include in output fasta file(or exclude if -n=True)
	-o, `-`-output_fasta_fp
		The output fasta file
	
	**[OPTIONAL]**
		
	-n, `-`-negate
		Negate the sample ID list (i.e., output sample ids not passed via -s) [default: False]


**Output:**

The script produces a fasta file containing containing only the specified SampleIDs.


**Examples:**

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE associated with sample ids S2, S3, S4 (-s). As always, sample IDs are case-sensitive:

::

	extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences  THAT ARE NOT (-n) associated with sample ids S2, S3, S4 (-s). As always, sample IDs are case-sensitive:

::

	extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4 -n


