.. _processing_illumina_data:

==========================
Processing Illumina Data
==========================

This document describes how to process Illumina sequencing data with QIIME. It covers a standard workflow beginning with fastq files, and take users through one pipeline for generation of an OTU table. Due to the huge amount of data in Illumina files, it's best to run the OTU picking through OTU table steps on a cluster. For a full HiSeq2000 run, this process can take up to 500 CPU hours.

This example illustrates how to use `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ to parse your Illumina output into a format that can be used by QIIME. 

Input file formats
^^^^^^^^^^^^^^^^^^
The preferred format in QIIME for Illumina data is fastq. This format is illustrated :ref:`here <fastq_format>`. Two other file formats can be used, but require additional steps detailed :ref:`here <other_file_formats>`. These steps convert other common formats to fastq, allowing you to begin processing your data with ``split_libraries_fastq.py``. 

To get information on the full set of parameters that can be passed to ``split_libraries_fastq.py``, review the `script documentation here <../scripts/split_libraries_fastq.html>`_.

Demultiplexing Illumina fastq with QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Demultiplexing and quality filtering one lane of Illumina fastq reads with QIIME
--------------------------------------------------------------------------------

To demultiplex your reads you'll need to know which of your fastq files corresponds to your barcodes, which corresponds to your amplicons, and you'll also need a metadata mapping file (format described `here <../documentation/file_formats.html#metadata-mapping-files>`_). To determine which are your barcode reads and which are your amplicon reads you can look at the first 4 lines of each file with the ``head`` command. For example::

	head -n 4 s_8_2_sequence.fastq
	
Will give you the first fastq record in s_8_2_sequence.fastq. If the length of the sequence (on the second line) corresponds to the length of your barcode, then this is your barcode read file. If not, check your other fastq file. Note: due to a sequencing artifact the barcode reads may be a single base longer than your actual barcodes. 

To demultiplex and quality filter (details on the quality filtering :ref:`here <illumina_quality_filtering>`) your fastq data, run the following command::

	split_libraries_fastq.py -i s_8_2_sequence.fastq -o sl_out/ -b s_8_1_sequence.fastq -m s_8_map.txt
	
In this example ``s_8_2_sequence.fastq`` contains my amplicon reads and ``s_8_1_sequence.fastq`` contains my barcode reads. My metadata mapping file is ``s_8_map.txt``, and the output will be written to ``./sl_out/``. This command can take tens of minutes to a couple of hours to run, depending on your computer. 

When the command completes you'll have three files in ``sl_out``:

 * seqs.fna : Your demultiplexed sequences.
 * split_libraries_log.txt : contains details on how many samples were assigned to each barcode, as well as details about sequences that were quality filter. It's a good idea to review this file after each run to confirm that the results look as expected.
 * histograms.txt : A histogram illustrating the distribution of sequence lengths after quality filtering/truncation.

After reviewing your split_libraries_log.txt you may choose to rerun with adjusted quality filter parameters, or continue to OTU picking. To continue to OTU picking move on to the next section.


Demultiplexing and quality filtering multiple lanes of Illumina fastq reads with QIIME
--------------------------------------------------------------------------------------

If you have multiple lanes of Illumina data that you want to demultiplex together, you'll use a similar command as for a single lane, but will specify per-lane amplicon read, barcode read, and mapping files. To do that, call the same command with comma-separated filepaths. Note that the order of the files must correspond for each of these parameters. For example, to demultiplex data from lanes 6, 7, and 8, your command would be::

	split_libraries_fastq.py -i s_6_2_sequence.fastq,s_7_2_sequence.fastq,s_8_2_sequence.fastq -o sl_out/ -b s_6_1_sequence.fastq,s_7_1_sequence.fastq,s_8_1_sequence.fastq -m s_6_map.txt,s_7_map.txts_8_map.txt
	
When the command completes you'll have the same output as when running on a single lane, but note that the data in the log and histogram files are broken down by lane. All of the sequence data will be in ``seqs.fna``.

Again, you may wish to review ``split_libraries_log.txt`` and adjust quality filtering parameters and rerun. When you're satisfied, you're read to move on to downstream analysis.

Using Illumina data in downstream analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have processed your Illumina data with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_, you can use the resulting files in all downstream QIIME scripts, including the workflow scripts.

A common way that we have processed Illumina data is using a closed reference OTU picking protocol, where OTUs are picked against a reference data set, and any reads that do not match a reference sequence at greater than or equal to 97% sequence identity are discarded. This has the benefit of providing a very strict quality filter (while studies to better understand quality filtering of Illumina amplicon data are in progress), but has the draw back of throwing away real sequences that do not match to know sequences.

To run a closed-reference OTU picking process on your demultiplexed Illumina data, you can use the `pick_reference_otus_through_otu_table.py <../scripts/pick_reference_otus_through_otu_table.html>`_ workflow. You can use your own reference database, or our greengenes OTUs reference set. To grab the latest version of these, click the "Most recent Greengenes OTUs" on the top-right corner of the `QIIME blog homepage <http://blog.qiime.org>`_. As of this writing the latest version is 4feb2011, but this is subject to change (hence not linking it directly from this page). Download this file and unzip it. Following from the steps above you can then pick reference otus on your demutliplexed data with the following command::

	pick_reference_otus_through_otu_table.py -i sl_out/seqs.fna -o ./ucrC/ -r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -t gg_otus_4feb2011/taxonomies/greengenes_tax.txt
	
In this example I am picking OTUs on ``sl_out/seqs.fna`` against the greengenes 97% reference OTU collection. The output is being written to ``./ucrC/`` (which stands for uclust_ref with the -C parameter, meaning 'closed reference'). Taxonomy strings from the Greengenes taxonomy is associated with each read based on their best match (determined by uclust_ref) in the reference set. 

This step will generate an OTU table, which is the input for a lot of the analyses possible with QIIME. For example, to generate interactive 3D UniFrac PCoA plots, you would run the command::

	beta_diversity_through_plots.py -i ucrC/uclust_ref_picked_otus/otu_table.biom -o bdiv/ -t gg_otus_4feb2011/trees/gg_97_otus_4feb2011.tre -m ./s_8_map.txt
	
Note that because we picked OTUs against a reference set, we can use the reference set phylogenetic tree for the UniFrac analysis. That is passed with ``-t`` in this example. To visualize the 3D UniFrac PCoA plots, you can open the ``bdiv/unweighted_unifrac_3d_continuous/unweighted_unifrac_pc_3D_PCoA_plots.html`` file that is generated in this analysis. This will launch the KiNG applet, and your 3D plots. These may take a little while to load depending on the quantity of data you have. (Improving these visualizations is something we're currently working on.)

.. _other_file_formats:

Processing non-fastq Illumina data with QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
QIIME supports several formats of non-fastq data, but the strategy is to convert from these formats to fastq. For that reason your analyses will be more convenient if you can get your sequencing center to provide data in fastq format (as supported by the Illumina CASAVA software).


Processing qseq files with QIIME
--------------------------------

You can convert qseq files to fastq files using the `process_qseq.py <../scripts/process_qseq.html>`_ script. 

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

You can convert iseq files to fastq files using the `process_iseq.py <../scripts/process_iseq.html>`_ script. Determine which of the following file types you have, and call the corresponding command.

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
	
Note that in the second example there are actually seven bases in the index field. If only six correspond to your barcode (and the remaining bases in e.g. a technical artifact) you can specify --barcode_length 6 (as done here) to extract only the first six bases of the barcode.

Once these steps are complete you'll have fastq files that can be passed to split_libraries_fastq.py.

Other topics
^^^^^^^^^^^^

.. _fastq_format:

Example fastq format (QIIME default)
------------------------------------

Example of amplicon read fastq::

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

Example of corresponding barcode read fastq::

	@M10_68:1:1:28680:29475#0/2
	ACTCACGGTATT
	+
	\_J\Sa^Y[ZYK
	@M10_68:1:1:19607:29475#0/2
	AGACTGAGTACT
	+
	PP\JJ\JQ`\RK
	@M10_68:1:1:22962:29475#0/2
	AGACGTGCAATT
	+
	^_aecceeeQ`[

.. _illumina_quality_filtering:

Quality filtering of Illumina data with QIIME
---------------------------------------------
A sequence is discarded if any of the following conditions are met:
	
	* The sequence contains one or more ``N`` bases, corresponding to ambiguous base calls (adjustable with the -n parameter).
	* The high-quality region of the sequence is less than 75 bases long (adjustable with the ``-p`` parameter), where high-quality regions is defined a stretch of bases containing no more than 1 (adjustable with the ``-r`` parameter) quality character less than ``B`` (i.e., any of ``@``, ``A``, and ``B`` are considered to be low quality scores, adjustable with the ``-q`` parameter). In other words, with the default parameter settings, the read is truncated at the base preceding the first low quality stretch, and the truncated sequence must be greater than or equal to 75 bases long to be retained. 
	* If barcode error correction is disabled and the barcode is not an exact match to a barcode in the mapping file (to disable this, pass ``-u``, which will cause the resulting sequences to be store with sample ID ``Unassigned``.)
	* If barcode error correction is enabled and the barcode is not correctable or within ``max_barcode_errors`` (specified on the command line with ``--max_barcode_errors``) of a good barcode. 


Processing paired-end read data with QIIME
------------------------------------------
QIIME can be used to parse single-end or paired-end read data from the Illumina platform. The downstream support for analysis of paired-end read data is currently limited. The parsed output is in standard fasta format, so all scripts (such as align_seqs.py and assign_taxonomy.py) can read it. However because there may be a 'big gap' between the 5' and 3' reads if the primers are distant in the sequence, or conversely because the reads may overlap if the primers are close, the reads are written to separate fasta files in separate runs. It is up to the user to merge these into a single file depending on how they wish to process the data (e.g., assemble over-lapping reads into contigs). 

If specific use cases become popular we will likely add support for them in QIIME. If you're interested in getting a specific workflow implemented you can contact us on the `QIIME Forum <http://forum.qiime.org>`_. The information we'll be interested in is an explanation of the workflow, and evidence that using paired-end reads improves results over using single-end reads alone. We are happy to share raw paired-end read Illumina data to facilitate such analyses.

To analyze demultiplex paired-end read data, run the split_libraries_fastq.py script on each read file independently. You can use the ``--rev_comp`` option on the reverse (3 prime) reads to reverse complement the reads so they'll be in the same orientation as the forward (5 prime) reads, if this is desired.

Barcode decoding
----------------
QIIME currently supports length 12 golay error-correcting barcodes as the default barcode type for ``split_libraries_fastq.py``. Non-golay codes of length 12, or barcodes of other lengths can be used by passing the length of the barcode with the ``--barcode_type`` option. No error recovery will be attempted in this case.

Barcode decoding slows down demultiplexing by a factor of approximately two (determined on a single HiSeq lane), but will be affected by how many barcodes actually need to be decoded. On this same HiSeq lane, barcode decoding recovered approximately 700,000 sequences (of approximately 60,000,000 input sequences) which had erroneous, correctable barcodes.



