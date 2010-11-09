.. _split_libraries_illumina:

.. index:: split_libraries_illumina.py

*split_libraries_illumina.py* -- Script for processing raw Illumina Genome Analyzer II data.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Script for parsing, library splitting, and quality filtering of raw Illumina Genome Analyzer II data.


**Usage:** :file:`split_libraries_illumina.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		The mapping filepath
	
	**[OPTIONAL]**
		
	-5, `-`-five_prime_read_fp
		The 5' read filepath [default: None]
	-3, `-`-three_prime_read_fp
		The 3' read filepath [default: None]
	-o, `-`-output_dir
		Output directory [default: ./]
	-u, `-`-store_unassigned
		Store seqs which can't be assigned to samples because of unknown barcodes [default: False]
	-q, `-`-quality_threshold
		Max base call error probability to consider high-quality (probability of base call being error, so values closer to 1 mean that the base call is more likely to be erroneous) [default: 1e-05]
	-r, `-`-max_bad_run_length
		Max number of consecutive low quality base calls allowed before truncating a read [default: 1; the read is trucated at thesecond low quality call]
	-p, `-`-min_per_read_length
		Min number of consecutive high quality base calls to includea read (per single end read) [default: 75]
	-n, `-`-sequence_max_n
		Maximum number of N characters allowed in a sequence to retain it -- this is applied after quality trimming, and is total over combined paired end reads if applicable [default: 0]
	-s, `-`-start_seq_id
		Start seq_ids as ascending integers beginning with start_seq_id[default: 0]
	`-`-rev_comp_barcode
		Reverse compliment barcodes before lookup[default: False]
	`-`-barcode_in_header
		Barcode is in header line (rather than beginning of sequence)[default: False]


**Output:**




**Parse paired-end read data (-5 and -3 provided), write output to s_1_seqs.fasta:**

::

	split_libraries_illumina.py -5 s_1_1_sequences.fasta -3 s_1_2_sequences.fasta -b barcode_map_6bp.txt

**Parse 5' read only (-5 only provided), write output to s_1_5prime_seqs.fasta:**

::

	split_libraries_illumina.py -5 s_1_1_sequences.fasta -b barcode_map_6bp.txt

**Parse 3' read only (-3 only provided), write output to  s_1_3prime_seqs.fasta:**

::

	split_libraries_illumina.py -3 s_1_2_sequences.fasta -b barcode_map_6bp.txt

**Parse multiple 5' read only files (multiple -5 values provided), write output to s_1_5prime_seqs.fasta, s_2_5primer_seqs.fasta:**

::

	split_libraries_illumina.py -5 s_1_1_sequences.fasta,s_2_1_sequences.fasta -b barcode_map_6bp.txt


