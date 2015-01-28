.. _processing_illumina_data:

========================================================================
Preparing raw Illumina data in different formats for use with QIIME
========================================================================

The `Illumina Overview Tutorial <./illumina_overview_tutorial.html>`_ describes how to work with raw Illumina sequence data with QIIME. That tutorial covers the case where you are starting with three specific input files: your `sample metadata mapping file  <../documentation/file_formats.html#metadata-mapping-files>`_ which contains the per-sample barcode sequences, a ``fastq`` file containing your amplicon sequence reads, and a corresponding ``fastq`` file containing the barcode reads for each amplicon sequence. If you don't have your data in this format, this document will help you get your data into this format for use with QIIME.

Demultiplexing and quality filtering of Illumina data is performed using `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_. The other scripts described in this document will help you to get your data into `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_.

.. note:: The preferred format in QIIME for Illumina data is ``fastq``. A detailed understanding of the ``fastq`` format is not necessary to use QIIME, but you can find details `here <http://scikit-bio.org/docs/latest/generated/skbio.io.fastq.html>`_ if you're interested. If your data is not in ``fastq`` format, this document can help you convert from some other less common formats to ``fastq``.

.. _illumina_quality_filtering:

Quality filtering of Illumina data with QIIME
=============================================

Quality filtering of Illumina data in QIIME is performed as part of `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_. Details on the quality filtering strategy in QIIME, including motivations for the default parameter settings, are described in `Bokulich et al. (2013) <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531572/>`_.

If you only wish to perform quality filtering of your fastq sequences (and not demultiplexing), you can use `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ with the ``--barcode_type 'not-barcoded'`` option.

Demultiplexing and quality filtering fastq Illumina data with QIIME
===================================================================

The standard approach for demultiplexing a single lane of Illumina fastq data with QIIME is illustrated in the `Illumina Overview Tutorial <./illumina_overview_tutorial.html>`_. Special cases are discussed here.

Starting with gzipped fastq files
---------------------------------

`split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ can work with either gzip-compressed (e.g., ``.fastq.gz``) or uncompressed (e.g. ``.fastq``) ``fastq`` files. You can pass compressed files to `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ in the same way that you pass uncompressed files.

Processing multiple lanes of Illumina fastq data with QIIME
-----------------------------------------------------------

If you have multiple lanes of Illumina data that you want to demultiplex together, you'll use a similar command as for a single lane (illustrated in the `Illumina Overview Tutorial <./illumina_overview_tutorial.html>`_), but will specify per-lane amplicon read, barcode read, and mapping files. To do that, call the same command with comma-separated filepaths. Note that the order of the files must correspond for each of these parameters. For example, to demultiplex data from lanes 6, 7, and 8, your command would be::

	split_libraries_fastq.py -i s_6_sequence.fastq,s_7_sequence.fastq,s_8_sequence.fastq -o sl_out/ -b s_6_barcode.fastq,s_7_barcode.fastq,s_8_barcode.fastq -m s_6_map.txt,s_7_map.txt,s_8_map.txt

When the command completes you'll have the same output as when running on a single lane, but note that the data in the log and histogram files are broken down by lane. All of the sequence data will be in ``seqs.fna``.

You may wish to review ``split_libraries_log.txt`` and adjust quality filtering parameters and rerun. When you're satisfied, you're read to move on to OTU picking, which is discussed in the `Illumina Overview Tutorial <./illumina_overview_tutorial.html>`_ and in `OTU picking strategies in QIIME <otu_picking.html>`_.

.. _other_file_formats:

Demultiplexing and quality filtering non-fastq Illumina data with QIIME
=======================================================================

Processing Illumina qseq files
------------------------------

You can convert qseq files to fastq files using the `process_qseq.py <../scripts/process_qseq.html>`_ script. The resulting fastq files can then be processed with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_.

Example qseq file for amplicon read::

	M10	68	1	1	28680	29475	0	1	AACGAAAGGCAGTTTTGGAAGTAGGCGAATTAGGGTAACGCATATAGGATGCTAATACAACGTGAATGAAGTACTGCATCTATGTCACCAGCTTATTACAGCAGCTTGTCATACATGGCCGTACAGGAAACACACATCATAGCATCACACGA	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	0
	M10	68	1	1	19607	29475	0	1	GACATAAGGGTGGTTAGTATACCGGCAAGGACGGGGTTACTAGTGACGTCCTTCCCCGTATGCCGGGCAATAATGTTTATGTTGGTTTCATGGTTTGGTCTAACTTTACCGCTACTAAATGCTGCGGATTGGTTTCGCTGAATCAGATTATT	Z__c\JQ`cc[[_[bfff[[`Qbdge_YYOOHO^cF[FUb_VHMHV`T`dBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	1
	M10	68	1	1	22962	29475	0	1	TAATCGAGCTCAACGCCCTGCATACGAAAAGACAGAATCTCTTGCAAGATGTTGGTGCGGTTAGCCAGCTGCTTATGGAAGCCAAGCATTGGGGATTGAGAAAGAGTAGAAATGCCACAAGCCTCAATAGCAGGTTTAAGAGCCTCGATACG	JJY````JO[`bab`b`bbaaaaa`\`a`OVT``]]`aa^aI\HMMMWWHHNNNGLL\`________\Z^]]^^^^^^GX]\QTXXZ[YZ^^XZ[Z^\Z^GW\^^\\^^^VZ\Y^^^^\\\\[^[\\\^VWYWWXWWZYZW^[X^\\Z^[TQ	0

Example qseq file for barcode read::

	M10	68	1	1	28680	29475	0	2	ACTCACGGTATTA	\_J\Sa^Y[ZYK`	0
	M10	68	1	1	19607	29475	0	2	AGACTGAGTACTA	PP\JJ\JQ`\RK^	1
	M10	68	1	1	22962	29475	0	2	AGACGTGCAATTA	^_aecceeeQ`[b	0

You'll need to know which of your reads files correspond to your barcodes and which correspond to your amplicons. To determine this you can look at the first line of representative files with the ``head`` command. For example::

	head -n 1 s_1_1_0001_qseq.txt

Will give you the first line s_1_1_0001_qseq.txt. If the length of the sequence (the 9th field) corresponds to the length of your barcode, then this is your barcode read file. If not, check a qseq file corresponding to another read number (e.g., s_1_2_0001_qseq.txt). Note that due to technical artifacts you may sometimes have a single extra base here, so for a length 12 barcode your sequence may be length 13.

You'll typically start here with a directory containing many qseq files. The process_qseq.py script therefore works on a directory, rather than a set of input files. In my example, the read 1 files correspond to my sequence reads and the read 2 files correspond to my barcode reads. To generate a single fastq file for the sequence reads from the qseq files, you can run the command::

	process_qseq.py -i ./ -o ./fastq/ -r 1

This specifies that the qseq files are in the current directory (``-i``), and the fastq should be written to ``./fastq/``. The ``-r 1`` specifies that I want to process the read one files (i.e., my amplicon reads).

To generate the barcode read fastq file you can run the following command::

	process_qseq.py -i ./ -o ./fastq/ -r 2 -b 12

This again specifies that the qseq files are in the current directory (``-i``), and the fastq should be written to ``./fastq/``. The ``-r 2`` specifies that I want to process the read two files (i.e., my barcode reads), and the ``-b 12`` specifies that I only want to extract the first twelve bases of these reads.

Once these steps are complete you'll have fastq files that can be passed to split_libraries_fastq.py.

Processing iseq files with QIIME
--------------------------------

You can convert iseq files to fastq files using the `process_iseq.py <../scripts/process_iseq.html>`_ script. The resulting fastq files can then be processed with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_. Determine which of the following file types you have, and call the corresponding command.

Example iseq with barcode in sequence (more common)::

	HWI-ST753_50:6:1101:15435:9071#0/1:ACCAGACGATGCTACGGAGGGAGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGCTTTGTAAGTTAGAGGTGAAAGCCTGGAGCTCAAC:gggggggfggdegggggggggggggggggggegggggggggegggggggeggcccccFUZSU_]]^^ggggggdggdgeeeccYacadcbeddceegggeeg
	HWI-ST753_50:6:1101:15446:9128#0/1:AGCTTAACAGCTTACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGCAGACGGAGAGGCAAGTCAGCTGTGAAAACTCCAGGCTTAAC:BBBBBBBBBBBB`_```_I^HM^`__`____I^^_`_`N``_______`__`___`_\_`G_^L^^^FDJTI^^^ZW^G^BBBBBBBBBBBBBBBBBBBBBB
	HWI-ST753_50:6:1101:15300:9134#0/1:ACCAGACGATGCTACGTAGGGGGCAAGCGTTGTCCGGATTTATTGGGCGTAAAGCGCGTGTAGGCGGCCAGGTAGGTCCGTTGTGAAAACTGGAGGCTTAAC:gggggggggcgcggggegggggeggfgggggggggggggggfggggggggggffMffa^cbbgggggggeggdedfb`dfeee`db^fffffge\geggdfg

To generate fastq from iseq files with tweleve base barcodes contained as the first bases of the sequence, call the following command::

	process_iseq.py -i s_6_1_sequences.txt,s_7_1_sequences.txt -o ./fastq/ -b 12


Example iseq with barcode in header (less common)::

	HWI-6X_9267:1:1:12:410#ACAGCTA/1:TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCATTTTAAGCCAGACGTGAAATCCCCGGGCTTAACCTGGGAACTG:abbb`aaa`^aa```ba`aaaabaaaabaaaa^[Y]^__a`abb`aaaa]Y\\_a[Y_a`a```a__]aaXT\`^\_]`a^^WSZ\JNY]^a`ORO^^`Y
	HWI-6X_9267:1:1:12:1762#ACATGAT/1:GACGGAGGATGCAAGTGTTATCCGGAATCACTGGGCGTAAAGCGTCTGTAGGTTGTTTGATAAGTCAACTGTTAAATCTTGAAGCTCAACTTCAAAATCG:aaaaaaaaabaaaaa_aaaaaa`aaaaaaaa`aa``a]aa```a^a^`\```\a`^aaa_\__]]_a_``^``a^^a^b[`SJN]Y_ZZ]^W___`_^U[
	HWI-6X_9267:1:1:12:1872#ACAGTTG/1:TACGGAGGGGGTTAGCGTTGTTCCGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGAAAGTTGGGGGGGAAATCCCGGGGCTCAACCCCGGACGTG:aaaaa_aaaa`[a_a`aaaa]a[MY``a\a`aaaaa_\]_\__[_]W]^[[U]aXRZ\W[J\KVTX]\YZZDVY]SUBBBBBBBBBBBBBBBBBBBBBBB

To generate fastq from iseq files with six base barcodes contained in the index field of the header, call the following command::

	process_iseq.py -i s_6_1_sequences.txt,s_7_1_sequences.txt -o ./fastq/ --barcode_length 6 --barcode_in_header

Note that in the second example there are actually seven bases in the index field. If only six correspond to your barcode (and the remaining bases are technical artifact, for example) you can specify --barcode_length 6 (as done here) to extract only the first six bases of the barcode.

Once these steps are complete you'll have fastq files that can be passed to `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_.

Processing paired-end read data with QIIME
==========================================

QIIME can be used to process single-end or paired-end read data from the Illumina platform. The primary script for merging paired-end read data in QIIME is `join_paired_ends.py <../scripts/join_paired_ends.html>`_. See the script documentation for more details. This is typically applied as a pre-processing step before running `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_.


Processing LEA-Seq data with QIIME
==================================

QIIME provides **beta support** for demultiplexing of sequences using the LEA-Seq protocol described in `Faith et al. (2013) <http://www.sciencemag.org/content/341/6141/1237439>`_. This can be performed using `split_libraries_lea_seq.py <../scripts/split_libraries_lea_seq.html>`_.


Working with already-demultiplexed fastq files
==============================================

Sequencing centers will often perform demultiplexing of sequences for you, delivering you one ``fastq`` file per sample. Several scripts are provided that facilitate using these files with QIIME:

 * `multiple_split_libraries_fastq.py <../scripts/multiple_split_libraries_fastq.html>`_: supports quality filtering of already-demultiplexed reads provided in multiple ``fastq`` files.
 * `multiple_join_paired_ends.py <../scripts/multiple_join_paired_ends.html>`_: supports joining paired ends for already-demultiplexed reads provided in multiple ``fastq`` files.
 * `multiple_extract_barcodes.py <../scripts/multiple_extract_barcodes.html>`_: supports removing barcode sequences from already-demultiplexed reads provided in multiple ``fastq`` files. This is useful when sequence barcodes were not removed as part of demultiplexing.
