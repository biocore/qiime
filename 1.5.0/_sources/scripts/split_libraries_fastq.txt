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
		
	`-`-store_qual_scores
		Store qual strings in .qual files [default: False]
	`-`-store_demultiplexed_fastq
		Write demultiplexed fastq files [default: False]
	`-`-retain_unassigned_reads
		Retain sequences which don't map to a barcode in the mapping file (sample ID will be "Unassigned") [default: False]
	-r, `-`-max_bad_run_length
		Max number of consecutive low quality base calls allowed before truncating a read [default: 3]
	-p, `-`-min_per_read_length_fraction
		Min number of consecutive high quality base calls to include a read (per single end read) as a fraction of the input read length [default: 0.75]
	-n, `-`-sequence_max_n
		Maximum number of N characters allowed in a sequence to retain it -- this is applied after quality trimming, and is total over combined paired end reads if applicable [default: 0]
	-s, `-`-start_seq_id
		Start seq_ids as ascending integers beginning with start_seq_id [default: 0]
	`-`-rev_comp_barcode
		Reverse compliment barcode reads before lookup [default: False]
	`-`-rev_comp_mapping_barcodes
		Reverse compliment barcode in mapping before lookup (useful if barcodes in mapping file are reverse compliments of golay codes)[default: False]
	`-`-rev_comp
		Reverse compliment sequence before writing to output file (useful for reverse-orientation reads) [default: False]
	-q, `-`-phred_quality_threshold
		The minimum acceptable Phred quality score (e.g., for Q20 and better, specify -q 20) [default: 3]
	`-`-last_bad_quality_char
		DEPRECATED: use -q instead. This method of setting is not robust to different versions of CASAVA.
	`-`-barcode_type
		The type of barcode used. This can be an integer, e.g. for length 6 barcodes, or golay_12 for golay error-correcting barcodes. Error correction will only be applied for golay_12 barcodes. [default: golay_12]
	`-`-max_barcode_errors
		Maximum number of errors in barcode [default: 1.5]


**Output:**




**Demultiplex and quality filter (at Phred Q20) one lane of Illumina fastq data and write results to ./slout_q20.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt -q20

**Demultiplex and quality filter (at Phred Q20) one lane of Illumina fastq data and write results to ./slout_q20. Store trimmed quality scores in addition to sequence data.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt --store_qual_scores -q20

**Demultiplex and quality filter (at Phred Q20) two lanes of Illumina fastq data and write results to ./slout_q20.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz,lane2_read1.fastq.gz -b lane1_barcode.fastq.gz,lane2_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt,map.txt --store_qual_scores -q20


