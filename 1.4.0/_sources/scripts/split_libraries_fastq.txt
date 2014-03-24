.. _split_libraries_fastq:

.. index:: split_libraries_fastq.py

*split_libraries_fastq.py* -- This script performs demultiplexing of Fastq sequence data where barcodes and sequences are contained in two separate fastq files (common on Illumina runs).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`split_libraries_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-sequence_read_fps
		The sequence read fastq files (comma-separated if more than one)
	-b, `-`-barcode_read_fps
		The barcode read fastq files (comma-separated if more than one)
	-o, `-`-output_dir
		Directory to store output files
	-m, `-`-mapping_fps
		Metadata mapping files (comma-separated if more than one)
	
	**[OPTIONAL]**
		
	`-`-retain_unassigned_reads
		Retain sequences which don't map to a barcode in the  mapping file (sample ID will be "Unassigned") [default: False]
	-r, `-`-max_bad_run_length
		Max number of consecutive low quality base calls allowed before truncating a read [default: 1; the read is trucated at thesecond low quality call]
	-p, `-`-min_per_read_length
		Min number of consecutive high quality base calls to includea read (per single end read) [default: 75]
	-n, `-`-sequence_max_n
		Maximum number of N characters allowed in a sequence to retain it -- this is applied after quality trimming, and is total over combined paired end reads if applicable [default: 0]
	-s, `-`-start_seq_id
		Start seq_ids as ascending integers beginning with start_seq_id[default: 0]
	`-`-rev_comp_barcode
		Reverse compliment barcode reads before lookup [default: False]
	`-`-rev_comp_mapping_barcodes
		Reverse compliment barcode in mapping before lookup (useful if barcodes in mapping file are reverse compliments of golay codes) [default: False]
	`-`-rev_comp
		Reverse compliment sequence before writing to output file (useful for reverse-orientation reads) [default: False]
	`-`-last_bad_quality_char
		The last character to be considered low quality (i.e., these character and those before it will be considered low quality base calls) [default: B]
	`-`-barcode_type
		The type of barcode used. This can be an integer, e.g. for length 6 barcodes, or golay_12 for golay error-correcting barcodes. Error correction will only be applied for golay_12 barcodes. [default: golay_12]
	`-`-max_barcode_errors
		Maximum number of errors in barcode [default: 1.5]


**Output:**




**Demultiplex and quality filter two lanes of Illumina run results and write results to ./sl_out/.:**

::

	split_libraries_fastq.py -i s_7_2_sequence.txt,s_8_2_sequence.txt -b s_7_1_sequence.txt,s_8_1_sequence.txt -o ./sl_out/ -m lane7_map.txt,lane8_map.txt


