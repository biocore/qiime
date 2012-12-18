.. _convert_fastaqual_fastq:

.. index:: convert_fastaqual_fastq.py

*convert_fastaqual_fastq.py* -- From a FASTA file and a matching QUAL file, generates a FASTQ file. From FASTQ file generates FASTA file and  matching QUAL file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

From a FASTA file and a matching QUAL file, generates a FASTQ file. A minimal FASTQ file omits the redundant sequence label on the quality scores; the quality scores for a sequence are assumed to follow immediately after the sequence with which they are associated. The output FASTQ file will be generated in the specified output directory with the same name as the input FASTA file, suffixed with '.fastq'. A FASTQ file will be split into FASTA and QUAL files, and generated in the designated output directory.


**Usage:** :file:`convert_fastaqual_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fasta_file_path
		Input FASTA or FASTQ file.
	
	**[OPTIONAL]**
		
	-q, `-`-qual_file_path
		Required input QUAL file if converting to FASTQ.
	-o, `-`-output_dir
		Output directory. Will be created if does not exist. [default: .]
	-c, `-`-conversion_type
		Type of conversion: fastaqual_to_fastq or fastq_to_fastaqual [default: fastaqual_to_fastq]
	-a, `-`-ascii_increment
		The number to add (subtract if coverting from FASTQ) to the quality score to get the ASCII character (or numeric quality score). [default: 33]
	-F, `-`-full_fasta_headers
		Include full FASTA headers in output file(s) (as opposed to merely the sequence label). [default: False]
	-b, `-`-full_fastq
		Include identifiers on quality lines in the FASTQ file (those beginning with a "+"). Irrelevant when converting from FASTQ. [default=False]
	-m, `-`-multiple_output_files
		Create multiple FASTQ files, one for each sample, or create multiple matching FASTA/QUAL for each sample. [default=False]


**Output:**

Outputs a complete or minimal FASTQ file, which omits the redundant sequence label on the quality scores, or splits FASTQ file into matching FASTA/QUAL files.


**Example:**

Using the input files seqs.fna and seqs.qual, generate seqs.fastq in the fastq_files directory:

::

	convert_fastaqual_fastq.py -f seqs.fna -q seqs.qual -o fastq_files/

**Example:**

Using input seqs.fastq generate fasta and qual files in fastaqual directory:

::

	convert_fastaqual_fastq.py -c fastq_to_fastaqual -f seqs.fastq -o fastaqual


