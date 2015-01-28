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

Working with alternate barcoding strategies
===========================================

In some cases, sequence barcodes are not provided in a separate file, or a dual barcoding strategy may have been applied during sequencing. QIIME supports working with these data by converting these formats into one that is compatible with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_. The examples below illustrate how to use `extract_barcodes.py <../scripts/extract_barcodes.html>`_ to accomplish this.

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

An example command to parse out the barcodes and the reads (with barcodes removed) to the output directory parsed_barcodes follows::

	extract_barcodes.py -f in_seqs.fastq --bc1_len 10 -o parsed_barcodes/ --input_type barcode_single_end

In the output directory, there should be a ``barcodes.fastq`` file with matching labels to the above file that contains the first 10 bases of the reads and the corresponding quality scores, and a ``reads.fastq`` containing the remaining bases and quality scores.

If one wanted to reverse complement the barcode before writing it, this command could be used::

	extract_barcodes.py -f in_seqs.fastq --bc1_len 10 --rev_comp_bc1 -o parsed_barcodes/ --input_type barcode_single_end

Processing two fastq files that each start with part of a barcode
-----------------------------------------------------------------

In this example, there are two reads, reads1.fastq and reads2.fastq. reads1.fastq starts with 6 bases of a barcode, and reads2.fastq has the remaining 8 bases of this barcode, which need to be pulled from the sequences and combined in an output ``barcodes.fastq`` file::

	extract_barcodes.py --input_type barcode_paired_end -f reads1.fastq -r reads2.fastq --bc1_len 6 --bc2_len 8 -o parsed_barcodes/

To change the order of the barcodes written, simply change the values/files of the ``-f``, ``-r``, ``--bc1_len``, and ``--bc2_len`` parameters::

	extract_barcodes.py --input_type barcode_paired_end -f reads2.fastq -r reads1.fastq --bc1_len 8 --bc2_len 6 -o parsed_barcodes/

The barcodes can each be reverse complemented before writing via the ``--rev_comp_bc1`` and ``--rev_comp_bc2`` parameters.

In some sequencing reactions, the orientation of the reads is random, so in the reads1.fastq and reads2.fastq example above, there would be a mixture of amplicons in different orientations. One solution is to attempt to orient the reads by looking for the forward and reverse primers in the sequences, and selectively writing out reads that match the orientation of the forward primer to reads1 and reads that match the reverse primer to reads2. One must supply a metadata mapping file that contains the LinkerPrimerSequence and ReversePrimer fields (all primers should be written in 5'->3' orientation). See `metadata mapping file format <../documentation/file_formats.html#metadata-mapping-files>`_ for more information. An example command to attempt read orientation::

	extract_barcodes.py --input_type barcode_paired_end -m mapping_file.txt -a -f reads1.fastq -r reads2.fastq --bc1_len 6 --bc2_len 8 -o parsed_barcodes/

In addition to the normal ``barcodes.fastq``, ``reads1.fastq``, and ``reads2.fastq`` files, there will be additional files labeled as ``_not_oriented.fastq`` where the primers could not be found. These are all written and have the barcodes processed in the order that the files were supplied in.

With paired sequences, one could run into issues of labels not matching-this can be suppressed with the ``--disable_header_match`` option, but one should be careful about making sure the reads are actually matched up if this option is used.

Two index/barcode reads and two fastq reads
-------------------------------------------

This situation can be treated as a special case of paired-end reads. One could supply the index files (labeled as index1.fastq, index2.fastq) and use the ``--input_type barcode_paired_end``::

	extract_barcodes.py --input_type barcode_paired_end -f index1.fastq -r index2.fastq --bc1_len 6 --bc2_len 6 -o parsed_barcodes/

The output ``barcodes.fastq`` file would be used for downstream processing, and the reads1 and reads2 files could be ignored.

Processing a single stitched read that contains barcodes on each end
--------------------------------------------------------------------

In this case, a single stitched fastq file (reads.fastq) is present that contains 6 base pair barcodes on each end. We can write the barcodes and the stitched read (sans barcodes) to the output directory parsed_barcodes with this command::

	extract_barcodes.py --input_type barcode_paired_stitched -f reads.fastq --bc1_len 6 --bc2_len 6 -o parsed_barcodes/

To switch the order of the barcodes (i.e., the barcode at the end of the read is written first), the ``--switch_bc_order`` option could be enabled. The barcodes can individually, or both, be reverse complemented in the same manner as described above. Additionally, one can attempt to orient the sequences as described in the paired-end read example. In this case, the orientation process will reverse complement the stitched reads to try and match the orientation of the forward primers.

Barcodes in fastq labels
------------------------

In some cases, sequencing centers will put the barcode sequences in the fastq labels before handing it off. It is important to look for the character immediately in front of the barcode and the length of the barcode::

	@MCIC-SOLEXA_0051_FC:1:1:4065:1039#CGATGT/1
	NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	+MCIC-SOLEXA_0051_FC:1:1:4065:1039#CGATGT/1
	KPPPQWWWWWQQ________BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

In this case, the "#" character is before the barcode, and the barcodes are 6 base pairs in length. To parse this example file, called in_seqs.fastq, this example command could be used::

	extract_barcodes.py --input_type barcode_in_label --char_delineator "#" -f in_seqs.fastq --bc1_len 6 -o parsed_barcodes/

A second fastq file could be passed (``-r``) if one had paired files with barcodes in the labels, and the parameters for changing barcode lengths or reverse complementing barcodes all apply.

Processing Joint Genome Institute fastq files
=============================================

The Joint Genome Institute (JGI) delivers a single fastq file to their users, which is not compatible with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_. To fix this issue, first we can split this file into forward and reverse reads and then extract the barcode from the forward file.

The fastq created by JGI alternates the forward and reverse reads in a single file and has the barcode in the header of the file, for example::

	@MISEQ03:64:000000000-A2H3D:1:1101:14358:1530 1:N:0:TCCACAGGAGT
	TNCAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGAAACTGACACTGAGGGGCGAAAGCGGGGGGGGCAAACG
	+
	?#5<????DDDDDDDDEEEEFFHHHHHHHHHHHHHHDCCHHFGDEHEH>CCE5AEEHHHHHHHHHHHHHHHHHFFFFHHHHHHEEADEEEEEEEEEEEEEEEEEEEEEEE?BEEEEEEEEEEEAEEEE0?A:?EE)8;)0ACEEECECCECAACEE?>)8CCC?CCA8?88ACC*A*::A??:0?C?.?0:?8884>'.''..'0?8C?C**0:0::?ECEE?############################
	@MISEQ03:64:000000000-A2H3D:1:1101:14358:1530 2:N:0:TCCACAGGAGT
	ACGGACTACAAGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGCACCTCAGTGTCAGTATCAGTCCAGGTGGTCGCCTTCGCCACTGGTGTTCCTTCCTATATCTACGCATTTCACCGCTACACAGGAAATTCCACCACCCTCTACCATACTCTAGCTTGTCAGTTTTGAATGCAGTTCCCAGGTTGAGCCCGGGGATTTCACATCCAACTTAACAACCCACATACCCGCCTTTTCGCCCAGGTAATCC
	+
	?????@@BDDBDDDBBCFFCFFHIIIIIIIIFGHHHHEHHHIIIHHHHHFHIIHIGHHIDGGHHHHIIIIICEFHIHHCDEHHHHHHFHHCFHDF?FHHFHHHFFDFFFDEDDD..=DDDE@<BFEEFCFFCECE==CACFE?*0:*CCAA?:*:*:0*A?A80:???A?*00:**0*1*:C??C?A?01*0?);>>'.8::A?###############################################
	@MISEQ03:64:000000000-A2H3D:1:1101:14206:1564 1:N:0:TCCACAGGAGT
	TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAGTCAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTGCGTTTGAAACTACAAGGCTAGAGTGTAGCAGAGGGGGGTAGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGGGGAGGAATACCAATGGCGAAGGCAGCCCCCGGGGTTAACACTGACGCCAAGGCACGAAAGCGGGGGGGGCAAACG
	+
	?????BB?DDDDDD@DDCEEFFH>EEFHHHHHHGHHHCEEFFDC5EECCCCCCDECEHF;?DFDDFHDDBBDF?CFDCCFEA@@::;EEEEEEEECBA,BBE?:>AA?CA*:**0:??A:8*:*0*0**0*:?CE?DD'..0????:*:?*0?EC*'.)4.?A***00)'.00*0*08)8??8*0:CEE*0:082.4;**?AEAA?#############################################
	@MISEQ03:64:000000000-A2H3D:1:1101:14206:1564 2:N:0:TCCACAGGAGT
	ACGGACTACAGGGGTTTCTAATCCTGTTTGCTCCCCACGCTTTCGTGCATGAGCGTCAGTGTTAACCCAGGGGGCTGCCTTCGCCATTGGTATTCCTCCACATCTCTACGCATTTCACTGCTACACGTGGAATTCTACCCCCCTCTGCTACACTCTAGCCTTGTAGTTTCAAACGCAGTTCCCAGGTTGAGCCCGGGGCTTTCACATCTGCCTTACAAAACCGCCTGCGCACGCTTTACGCCCCGTAATTC
	+
	?????@@BDDBDDD?@CFFFFFHHHHHFFHHHHHHHHHHH@FFHEDFFH@FHBGCDHHHBFHHHHHHHEHHHHDCCEFFDFFFEE:=?FF?DFDFDFFF==BEE=DBDDEEEEEB,4??EE@EEE,3,3*3,?:?*0ACCEDD88:***?*0:*0***0*?C?00:AE:?EE:*A8'.?:CAA?A80*0*??AA88;28;C##################################################

To split the input file in forward and reverse we will use `extract_reads_from_interleaved_file.py <../scripts/extract_reads_from_interleaved_file.html>`_::

	extract_reads_from_interleaved_file.py -i reads_to_extract.fastq -o extracted_reads

Then the next step is to create the barcode (index) fastq file by using the forward read generated in the previous step and `extract_barcodes.py <../scripts/extract_barcodes.html>`_. Note how we are using the parameter ``-l 11`` because the primers in this example have that length, see header in the example::

	extract_barcodes.py -f extracted_reads/forward_reads.fastq -o barcode_reads -c barcode_in_label --char_delineator '1:N:0:' -l 11
