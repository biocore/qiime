.. _extracting_barcodes:

=====================================================================================
Extracting Barcodes from fastq data for compatibility with `split_libraries_fastq.py`
=====================================================================================

Introduction
------------

To use barcoded demultiplexing with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_, one must supply a fastq file with the barcodes, and another fastq file with the reads (the labels must match between these files). Please see `this <./processing_illumina_data.html#fastq-format>`_ for more details.

In some cases, alternative data formats are generated due to how primer constructs are designed, or reads are processed before using them in QIIME. If one's data matches the examples below, `extract_barcodes.py <../scripts/extract_barcodes.html>`_ can be used to parse out the barcodes and reads into a format compatible with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_.

Processing a single fastq file that starts with the barcode sequence
--------------------------------------------------------------------

In this example, there is a single fastq file, and the first 10 base pairs of the fastq file (in_seqs.fastq) are the barcode

Example of raw data for the file in_seqs.fastq::

	@M10_68:1:1:28680:29475#0/1
	AACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACGA
	+
	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
	@M10_68:1:1:19607:29475#0/1
	GACATAAGGGTGGTTAGTATACCGGCAAGGACGGGGTTACTAGTGACGTCCTTCCCCGTATGCCGGGCAATAATGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCTGCGGATTGGTTTCGCTGAATCAGATTATT
	+
	Z__c\JQ`cc[[_[bfff[[`Qbdge_YYOOHO^cF[FUb_VHMHV`T`dBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
	@M10_68:1:1:22962:29475#0/1
	TAATCGAGCTCAACGCCCTGCATACGAAAAGACAGAATCTCTTGCAAGATGTTGGTGCGGTTAGCCAGCTGCTTATGGAAGCCAAGCATTGGGGATTGAGAAAGAGTAGAAATGCCACAAGCCTCAATAGCAGGTTTAAGAGCCTCGATACG
	+
	JJY````JO[`bab`b`bbaaaaa`\`a`OVT``]]`aa^aI\HMMMWWHHNNNGLL\`________\Z^]]^^^^^^GX]\QTXXZ[YZ^^XZ[Z^\Z^GW\^^\\^^^VZ\Y^^^^\\\\[^[\\\^VWYWWXWWZYZW^[X^\\Z^[TQ

An example command to parse out the barcodes and the reads (with barcodes removed) to the output directory parsed_barcodes would be this::

	extract_barcodes.py -f in_seqs.fastq --bc1_len 10 -o parsed_barcodes/ --input_type barcode_single_end

In the output directory, there should be a ``barcodes.fastq`` file with matching labels to the above file that contains the first 10 bases of the reads and the corresponding quality scores, and a ``reads.fastq`` containing the remaining bases and quality scores.

If one wanted to reverse complement the barcode before writing it, this command could be used::

	extract_barcodes.py -f in_seqs.fastq --bc1_len 10 --rev_comp_bc1 -o parsed_barcodes/ --input_type barcode_single_end

Processing two fastq files that each start with part of a barcode
-----------------------------------------------------------------

In this example, there are two reads, reads1.fastq and reads2.fastq. reads1.fastq starts with 6 bases of a barcode, and reads2.fastq has the remaining 8 bases of this barcode, which need to be pulled from the sequences and combined in an output ``barcodes.fastq`` file::

	extract_barcodes.py --input_type barcode_paired_end -f reads1.fastq -F reads2.fastq --bc1_len 6 --bc2_len 8 -o parsed_barcodes/
	
To change the order of the barcodes written, simply change the values/files of the ``-f``, ``-F``, ``--bc1_len``, and ``--bc2_len`` parameters::

	extract_barcodes.py --input_type barcode_paired_end -f reads2.fastq -F reads1.fastq --bc1_len 8 --bc2_len 6 -o parsed_barcodes/

The barcodes can each be reverse complemented before writing via the ``--rev_comp_bc1`` and ``--rev_comp_bc2`` parameters.

In some sequencing reactions, the orientation of the reads is random, so in the reads1.fastq and reads2.fastq example above, there would be a mixture of amplicons in different orientations. One solution is to attempt to orient the reads by looking for the forward and reverse primers in the sequences, and selectively writing out reads that match the orientation of the forward primer to reads1 and reads that match the reverse primer to reads2. One must supply a metadata mapping file that contains the LinkerPrimerSequence and ReversePrimer fields (all primers should be written in 5'->3' orientation). See `metadata mapping file format <../documentation/file_formats.html#metadata-mapping-files>`_ for more information. An example command to attempt read orientation::

	extract_barcodes.py --input_type barcode_paired_end -m mapping_file.txt -a -f reads1.fastq -F reads2.fastq --bc1_len 6 --bc2_len 8 -o parsed_barcodes/

In addition to the normal ``barcodes.fastq``, ``reads1.fastq``, and ``reads2.fastq`` files, there will be additional files labeled as ``_not_oriented.fastq`` where the primers could not be found. These are all written and have the barcodes processed in the order that the files were supplied in.

With paired sequences, one could run into issues of labels not matching-this can be suppressed with the ``--disable_header_match`` option, but one should be careful about making sure the reads are actually matched up if this option is used.

Two index/barcode reads and two fastq reads
-------------------------------------------

This situation can be treated as a special case of paired-end reads. One could supply the index files (labeled as index1.fastq, index2.fastq) and use the ``--input_type barcode_paired_end``::

	extract_barcodes.py --input_type barcode_paired_end -f index1.fastq -F index2.fastq --bc1_len 6 --bc2_len 6 -o parsed_barcodes/

The output ``barcodes.fastq`` file would be used for downstream processing, and the reads1 and reads2 files could be ignored.

Processing a single stitched read that contains barcodes on each end
--------------------------------------------------------------------

In this case, a single stitched fastq file (reads.fastq) is present that contains 6 base pair barcodes on each end. We can write the barcodes and the stitched read (sans barcodes) to the output directory parsed_barcodes with this command::

	extract_barcodes.py --input_type barcode_paired_stitched -f reads.fastq --bc1_len 6 --bc2_len 6 -o parsed_barcodes/

If one wanted the order of the barcodes switched (the barcode at the end of the read written first), the ``--switch_bc_order`` option could be enabled. The barcodes can individually, or both, be reverse complemented in the same manner as described above. Additionally, one can attempt to orient the sequences as described in the paired-end read example. In this case, the orientation process will reverse complement the stitched reads to try and match the orientation of the forward primers.

Barcodes in fastq labels
------------------------

In some cases, sequencing centers will put the barcode sequences in the fastq labels before handing it off. It is important to look for the character immediately in front of the barcode and the length of the barcode::

	@MCIC-SOLEXA_0051_FC:1:1:4065:1039#CGATGT/1
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	+MCIC-SOLEXA_0051_FC:1:1:4065:1039#CGATGT/1
	KPPPQWWWWWQQ________BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

In this case, the "#" character is before the barcode, and the barcodes are 6 base pairs in length. To parse this example file, called in_seqs.fastq, this example command could be used::

	extract_barcodes.py --input_type barcode_in_label --char_delineator "#" -f in_seqs.fastq --bc1_len 6 -o parsed_barcodes/
	
A second fastq file could be passed (``-F``) if one had paired files with barcodes in the labels, and the parameters for changing barcode lengths or reverse complementing barcodes all apply.

Notes for post-demultiplexing
-----------------------------

In many of these cases, the primer sequences (forward and potentially reverse) will remain in the sequences. It is standard practice to remove these from the sequences before clustering or other analyses. The `quality_filter_fastq.py` script can remove forward and optionally, reverse primers, but one should disable the other quality filtering settings (e.g. homopolymer checks, which are issues with the 454 platform).