.. _convert_fastaqual_to_fastq:

.. index:: convert_fastaqual_to_fastq.py

*convert_fastaqual_to_fastq.py* -- From a FASTA file and a matching QUAL file, generates a minimal FASTQ file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

From a FASTA file and a mathcing QUAL file, generates a minimal FASTQ file. A minimal FASTQ file omits the redundtant sequence label on the quality scores; the quality scores for a sequence are assumed to follow immediately after the sequence with which they are associated. The output FASTQ file will be generated in the specified output directory with the same name as the input FASTA file, suffixed with '.fastq'


**Usage:** :file:`convert_fastaqual_to_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fasta_fp
		Input FASTA file.
	-q, `-`-qual_fp
		Input QUAL file.
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Output directory. Will be created if does not exist. [default: .]
	-a, `-`-ascii_increment
		The number to add to the quality score to get the ASCII character. [default: 33]
	-F, `-`-full_fasta_headers
		Include full FASTA headers in FASTQ file (as opposed to merely the sequence label). [default: False]
	-b, `-`-full_fastq
		Include identifiers on quality lines in the FASTQ file (those beginning with a "+" [default=False]
	-m, `-`-multiple_output_files
		Create multiple FASTQ files, one for each sample. [default=False]


**Output:**

Outputs a minimal FASTQ file, which omits the redundant sequence label on the quality scores.


**Example:**

Using the input files seqs.fna and seqs.qual, generate seqs.fastq in the fastq_files directory:

::

	python fastaQualToFastq_script.py -f seqs.fna -q seqs.qual -o fastq_files/


