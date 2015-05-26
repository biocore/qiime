.. _subsample_fasta:

.. index:: subsample_fasta.py

*subsample_fasta.py* -- Randomly subsample sequences from a given fasta file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Subsample the seqs.fna file, randomly select 5% of the sequences:


**Usage:** :file:`subsample_fasta.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-p, `-`-percent_subsample
		Specify the percentage of sequences to subsample
	
	**[OPTIONAL]**
		
	-o, `-`-output_fp
		The output filepath


**Output:**




**Example:**

 

::

	 subsample_fasta.py -i seqs.fasta -p 0.05


