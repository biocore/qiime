.. _split_libraries_fastq:

.. index:: split_libraries_fastq.py

*split_libraries_fastq.py* -- This script performs demultiplexing of Fastq sequence data where barcodes and sequences are contained in two separate fastq files (common on Illumina runs).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The Illumina MiSeq platform commonly outputs sequence reads and barcode reads in separate files. For example:

.. note::
- "Example_S0_L001_R1_001.fastq.gz": 	R1 is read one, the forward read
- "Example_S0_L001_I1_001.fastq.gz": 	I1 is the index, the barcode read

When paired end reads are generated on an Illumina machine, they can be paired with `join_paired_ends.py <http://qiime.org/scripts/join_paired_ends.html>`_ and produce the following output:

.. note::
- "fastqjoin.join.fastq":		the result of pairing the forward and reverse reads
- "fastqjoin.join_barcodes.fastq": 	the updated barcodes for these newly paired reads

In both these cases, this script can be used to demultiplex the reads.fastq file using the barcodes.fastq file and `mapping file <http://qiime.org/documentation/file_formats.html#mapping-file-overview>`_.



**Usage:** :file:`split_libraries_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-sequence_read_fps
		The sequence read fastq files (comma-separated if more than one)
	-o, `-`-output_dir
		Directory to store output files
	
	**[OPTIONAL]**
		
	-m, `-`-mapping_fps
		Metadata mapping files (comma-separated if more than one) [default: None]
	-b, `-`-barcode_read_fps
		The barcode read fastq files (comma-separated if more than one) [default: None]
	`-`-store_qual_scores
		Store qual strings in .qual files [default: False]
	`-`-sample_ids
		Comma-separated list of samples ids to be applied to all sequences, must be one per input file path (used when data is not multiplexed) [default: None]
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
		Reverse complement barcode reads before lookup [default: False]
	`-`-rev_comp_mapping_barcodes
		Reverse complement barcode in mapping before lookup (useful if barcodes in mapping file are reverse complements of golay codes) [default: False]
	`-`-rev_comp
		Reverse complement sequence before writing to output file (useful for reverse-orientation reads) [default: False]
	-q, `-`-phred_quality_threshold
		The maximum unacceptable Phred quality score (e.g., for Q20 and better, specify -q 19) [default: 3]
	`-`-last_bad_quality_char
		DEPRECATED: use -q instead. This method of setting is not robust to different versions of CASAVA.
	`-`-barcode_type
		The type of barcode used. This can be an integer, e.g. for length 6 barcodes, or "golay_12" for golay error-correcting barcodes. Error correction will only be applied for "golay_12" barcodes. If data is not barcoded, pass "not-barcoded". [default: golay_12]
	`-`-max_barcode_errors
		Maximum number of errors in barcode [default: 1.5]
	`-`-phred_offset
		The ascii offset to use when decoding phred scores (either 33 or 64). Warning: in most cases you don't need to pass this value [default: determined automatically]


**Output:**




**Demultiplex and quality filter (at Phred >= Q20) one lane of Illumina fastq data and write results to ./slout_q20.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt -q 19

**Demultiplex and quality filter (at Phred >= Q20) one lane of Illumina fastq data and write results to ./slout_q20. Store trimmed quality scores in addition to sequence data.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt --store_qual_scores -q 19

**Demultiplex and quality filter (at Phred >= Q20) two lanes of Illumina fastq data and write results to ./slout_q20.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz,lane2_read1.fastq.gz -b lane1_barcode.fastq.gz,lane2_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ -m map.txt,map.txt --store_qual_scores -q 19

**Quality filter (at Phred >= Q20) one non-multiplexed lane of Illumina fastq data and write results to ./slout_single_sample_q20.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz --sample_ids my.sample.1 -o slout_single_sample_q20/ -q 19 --barcode_type 'not-barcoded'

**Quality filter (at Phred >= Q20) two non-multiplexed lanes of Illumina fastq data with different samples in each and write results to ./slout_not_multiplexed_q20.:**

::

	split_libraries_fastq.py -i lane1_read1.fastq.gz,lane2_read1.fastq.gz --sample_ids my.sample.1,my.sample.2 -o slout_not_multiplexed_q20/ -q 19 --barcode_type 'not-barcoded'


