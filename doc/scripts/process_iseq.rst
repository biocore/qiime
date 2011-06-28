.. _process_iseq:

.. index:: process_iseq.py

*process_iseq.py* -- Given a directory of per-swath qseq files, this script generates a single fastq per lane.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`process_iseq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		The input filepaths (either iseq or gzipped iseq format; comma-separated if more than one). See Processing Illumina Data tutorial for a description of the iseq file type.
	-o, `-`-output_dir
		The output directory
	-b, `-`-barcode_length
		Length of the barcode
	
	**[OPTIONAL]**
		
	`-`-barcode_in_header
		Pass if barcode is in the header index field (rather than at the beginning of the sequence)
	`-`-barcode_qual_c
		If no barcode quality string is available, score each base with this quality [default: b]


**Output:**




Generate fastq files from lanes 1 and 2 (read 1 data) where barcodes are contained as the first tweleve bases of the sequences.

::

	process_qseq.py -i ./s_1_1_sequence.txt,./s_2_1_sequence.txt -b 12 -o ./fastq/

Generate fastq files from the gzipped lanes 1 and 2 (read 1 data) where barcodes are contained as the first tweleve bases of the sequences.

::

	process_qseq.py -i ./s_1_1_sequence.txt.gz,./s_2_1_sequence.txt.gz -b 12 -o ./fastq/


