.. _extract_seqs_by_sample_id:

.. index:: extract_seqs_by_sample_id.py

*extract_seqs_by_sample_id.py* -- Extract sequences based on the SampleID
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script creates a fasta file which will contain only sequences that ARE associated with a set of sample IDs, OR all sequences that are NOT associated with a set of sample IDs (-n)


**Usage:** :file:`extract_seqs_by_sample_id.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-o, `-`-output_fasta_fp
		The output fasta file
	
	**[OPTIONAL]**
		
	-n, `-`-negate
		Negate the sample ID list (i.e., output sample ids not passed via -s) [default: False]
	-s, `-`-sample_ids
		Comma-separated sample_ids to include in output fasta file (or exclude if --negate), or string describing mapping file states defining sample ids (mapping_fp must be provided for the latter)
	-m, `-`-mapping_fp
		The mapping filepath


**Output:**

The script produces a fasta file containing containing only the specified SampleIDs.


**Examples:**

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE associated with sample ids S2, S3, S4 (-s). As always, sample IDs are case-sensitive:

::

	extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs_by_sample.fasta -s S2,S3,S4

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE NOT (-n) associated with sample ids S2, S3, S4 (-s). As always, sample IDs are case-sensitive:

::

	extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs_by_sample_negated.fasta -s S2,S3,S4 -n

Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE associated with sample ids whose "Treatment" value is "Fast" in the mapping file:

::

	extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs_by_mapping_field.fasta -m map.txt -s "Treatment:Fast" 


