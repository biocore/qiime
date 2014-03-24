.. _blast_wrapper:

.. index:: blast_wrapper.py

*blast_wrapper.py* -- Blast Interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script is a functionally-limited interface to the qiime.util.qiime_blast_seqs function, primarily useful for testing purposes. Once that function has been integrated into qiime as the primary blast interface it will move to PyCogent. An expanded version of this command line interface may replace the script functionality of cogent.app.blast at that point.


**Usage:** :file:`blast_wrapper.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-r, `-`-refseqs_fp
		Path to blast database as a fasta file
	
	**[OPTIONAL]**
		
	-n, `-`-num_seqs_per_blast_run
		Number of sequences passed to each blast call - useful for very large sequence collections [default: 1000]


**Output:**

This is a utility program, which returns BLAST results.


**Example:**

Blast all sequences in inseqs.fasta (-i) against a BLAST db constructed from refseqs.fasta (-r).

::

	blast_wrapper.py -i inseqs.fasta -r refseqs.fasta


