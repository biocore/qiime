.. _truncate_fasta_qual_files:

.. index:: truncate_fasta_qual_files.py

*truncate_fasta_qual_files.py* -- Generates filtered fasta and quality score files by truncating at the specified base position.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This module is designed to remove regions of poor quality in
454 sequence data.  Drops in quality can be visualized with the
`quality_scores_plot.py <./quality_scores_plot.html>`_ module.  The base position specified will
be used as an index to truncate the sequence and quality scores, and
all data at that base position and to the end of the sequence will be
removed in the output filtered files.


**Usage:** :file:`truncate_fasta_qual_files.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fasta_fp
		Input fasta filepath to be truncated.
	-q, `-`-qual_fp
		Input quality scores filepath to be truncated.
	-b, `-`-base_pos
		Nucleotide position to truncate the fasta and quality score files at.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Output directory.  Will be created if does not exist.  [default: .]


**Output:**

Filtered versions of the input fasta and qual file (based on input name with '_filtered' appended) will be generated in the output directory


**Example:**

Truncate the input fasta and quality files at base position 100, output to the filtered_seqs directory:

::

	truncate_fasta_qual_files.py -f seqs.fna -q seqs.qual -b 100 -o filtered_seqs/


