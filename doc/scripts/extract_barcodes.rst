.. _extract_barcodes:

.. index:: extract_barcodes.py

*extract_barcodes.py* -- This script is designed to format fastq sequence and barcode data so they are compatible with split_libraries_fastq.py (see http://qiime.org/tutorials/processing_illumina_data.html).
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

A variety of data formats are possible, depending upon how one utilized sequencing primers, designed primer constructs (e.g., partial barcodes on each end of the read), or processed the data (e.g., barcodes were put into the sequence labels rather than the reads). See various input examples below.


**Usage:** :file:`extract_barcodes.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-fastq1
		Input fastq filepath. This file is considered read 1.
	
	**[OPTIONAL]**
		
	-r, `-`-fastq2
		Input fastq filepath. This file is considered read 2. [default: None]
	-o, `-`-output_dir
		Directory prefix for output files [default: .]
	-c, `-`-input_type
		Specify the input type. barcode_single_end: Input is a single fastq file, that starts with the barcode sequence. barcode_paired_end: Input is a pair of fastq files (--fastq1 and --fastq2) that each begin with a barcode sequence. The barcode for fastq1 will be written first, followed by the barcode from fastq2. barcode_paired_stitched: Input is a single fastq file that has barcodes at the beginning and end. The barcode from the beginning of the read will be written first followed by the barcode from the end of the read, unless the order is switched with --switch_bc_order. barcode_in_label: Input is a one (--fastq1) or two (--fastq2) fastq files with the barcode written in the labels. [default: barcode_single_end]
	-l, `-`-bc1_len
		Specify the length, in base pairs, of barcode 1. This applies to the --fastq1 file and all options specified by --input_type [default: 6]
	-L, `-`-bc2_len
		Specify the length, in base pairs, of barcode 2. This applies to the --fastq2 file and options "barcode_paired_end", "barcode_paired_stitched", and "barcode_in_label" for the --input_type [default: 6]
	`-`-rev_comp_bc1
		Reverse complement barcode 1 before writing [default: False]
	`-`-rev_comp_bc2
		Reverse complement barcode 2 before writing [default: False]
	-s, `-`-char_delineator
		Character in fastq label that should immediately precede the barcode sequence. The length of the barcode is specified by the --bc1_len (and optionally --bc2_len if paired end files are used) parameter. [default: :]
	`-`-switch_bc_order
		Reverse barcode order written when using the -c barcode_paired_stitched option. [default: False]
	-m, `-`-mapping_fp
		Filepath of mapping file. NOTE: Must contain a header line indicating SampleID in the first column and BarcodeSequence in the second, LinkerPrimerSequence in the third and a ReversePrimer column before the final Description column. Needed for --attempt_read_orientation option. [default: None]
	-a, `-`-attempt_read_reorientation
		Will attempt to search for the forward and reverse primer in the read and adjust the sequence orientation to match the orientation of the forward primer. An exact match for the  forward and reverse complemented versions of the primers are tested for, and sequences are reverse complemented, if necessary, before writing. Sequences without an exact match are written to a separate output fastq file, labeled as _no_primer_match.fastq. [default: False]
	-d, `-`-disable_header_match
		Enable this option to suppress header matching between input fastq files.[default: False]


**Output:**

In the output directory, there will be fastq files (barcode file, and one or two reads files)


**Parse barcodes of 12 base pairs from the beginning of a single read. Will create an output fastq file of the barcodes and an output file of the reads supplied with the barcodes removed.:**

::

	extract_barcodes.py -f inseqs.fastq -c barcode_single_end --bc1_len 12 -o processed_seqs

**Parse barcodes of 12 base pairs from the beginning of a single read, reverse complement the barcodes before writing. Will create an output fastq file of the barcodes and an output file of the reads supplied with the barcodes removed:**

::

	extract_barcodes.py -f inseqs.fastq -c barcode_single_end --bc1_len 12 -o processed_seqs --rev_comp_bc1

**Parse barcodes of 6 base pairs from the beginning of paired reads. Will create an output fastq file of the barcodes and an output file of each of the reads supplied with the barcodes removed. The order of the barcodes written is determined by the order of the files passed (-f is written first, followed by -r):**

::

	extract_barcodes.py -f inseqs_R1.fastq -r inseqs_R2.fastq -c barcode_paired_end --bc1_len 6 --bc2_len 6 -o processed_seqs

**Parse barcodes of 6 base pairs from the beginning of paired reads, attempt to orient reads based upon detection of forward and reverse primers in the mapping file. Will create an output fastq file of the barcodes and an output file of each of the reads supplied with the barcodes removed. The order of the barcodes written is determined by the order of the files passed (-f is written first, followed by -r):**

::

	extract_barcodes.py -f inseqs_R1.fastq -r inseqs_R2.fastq -c barcode_paired_end --map_fp mapping_data.txt --attempt_read_orientation --bc1_len 6 --bc2_len 6 -o processed_seqs

**Parse barcodes of 6 base pairs from the beginning, 8 base pairs at the end of a stitched read. Will create an output fastq file of the barcodes and an output fastq file of the stitched read supplied with the barcodes removed. The barcode at the beginning of the stitched read is written first, followed by the barcode at the end, unless reversed by the --switch_bc_order option is used:**

::

	extract_barcodes.py -f inseqs_R1.fastq -c barcode_paired_stitched --bc1_len 6 --bc2_len 8 -o processed_seqs

**Parse barcodes of 12 base pairs from labels of the input fastq file. Example label (note that the desired character preceding the barcode is '#'): @MCIC-SOLEXA_0051_FC:1:1:14637:1026#CGATGTGATTTC/1 This will create an output fastq file of the barcodes (no other sequence are written). A second file with barcodes in the label can be passed with -r, and if this is done, the combined barcodes from -f and -r will be written together:**

::

	extract_barcodes.py -f inseqs_R1.fastq -c barcode_in_label --char_delineator '#' --bc1_len 12 -o processed_seqs


