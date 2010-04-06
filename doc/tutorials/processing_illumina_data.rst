.. _processing_illumina_data:

==========================
Processing Illumina Data
==========================

This document describes how to process Illumina sequencing data with QIIME. As the Illumina platform is just beginning to be used for community sequencing, the Illumina support in QIIME is still in development status. If you run into issues because your Illumina output is different from what we assume here, we recommend getting in touch via qiime.help@colorado.edu. Full support for analyzing community sequencing data generated on the Illumina platform is one of our goals, so we're willing to help out.

QIIME can be used to parse single-end or paired-end read data from the Illumina platform. The downstream support for analysis of paired-end read data is currently more limited. The parsed output is in standard fasta format, so all scripts (such as align_seqs.py and assign_taxonomy.py) can read it. However because there may be a 'big gap' between the 5' and 3' reads if the primers are distant in the sequence, or conversely because the reads may overlap if the primers are close, your mileage with the downstream tools may vary. 

This example illustrates how to used `split_libraries_illumina.py <../scripts/split_libraries_illumina.html>`_ to parse your Illumina output into a format that can be used by QIIME. 

Input file formats
-------------------------

The input to split_libraries_illumina.py is one or more ``sequence`` files, and a standard QIIME mapping file. For example, 10 lines extracted from a ``s_1_1_sequence.txt`` file (where the first ``1`` refers to the lane, and the second ``1`` refers to the read, indicating the 5' read in this case) are::


	HWI-6X_9267:1:1:12:410#ACAGCTA/1:TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGCATTTTAAGCCAGACGTGAAATCCCCGGGCTTAACCTGGGAACTG:abbb`aaa`^aa```ba`aaaabaaaabaaaa^[Y]^__a`abb`aaaa]Y\\_a[Y_a`a```a__]aaXT\`^\_]`a^^WSZ\JNY]^a`ORO^^`Y
	HWI-6X_9267:1:1:12:1762#ACATGAT/1:GACGGAGGATGCAAGTGTTATCCGGAATCACTGGGCGTAAAGCGTCTGTAGGTTGTTTGATAAGTCAACTGTTAAATCTTGAAGCTCAACTTCAAAATCG:aaaaaaaaabaaaaa_aaaaaa`aaaaaaaa`aa``a]aa```a^a^`\```\a`^aaa_\__]]_a_``^``a^^a^b[`SJN]Y_ZZ]^W___`_^U[
	HWI-6X_9267:1:1:12:1872#ACAGTTG/1:TACGGAGGGGGTTAGCGTTGTTCCGAATTACTGGGCGTAAAGCGCGCGTAGGCGGATTAGAAAGTTGGGGGGGAAATCCCGGGGCTCAACCCCGGACGTG:aaaaa_aaaa`[a_a`aaaa]a[MY``a\a`aaaaa_\]_\__[_]W]^[[U]aXRZ\W[J\KVTX]\YZZDVY]SUBBBBBBBBBBBBBBBBBBBBBBB


You'll notice that each line contains seven ``:`` delimited fields. These are:
	
	#. Machine name
	#. Channel/lane number
	#. Tile number
	#. X position
	#. Y position
	#. nucleotide sequence
	#. quality score

In this case, the ``Y position`` field additionally contains the reverse compliment of the barcode sequence as the first ``B`` bases following the ``#`` sign, where ``B`` refers to the barcode length as determined from the mapping file. If your data looks different from this, get in touch with qiime.help@colorado.edu, and we can try to work on custom solutions until standards are developed in this area.

One important thing to note is that our sequencing primers were developed to avoid sequencing the PCR primers. So unlike standard 454 data, there are no barcode, linker, or primer sequences that need to be extracted from these sequences. As a consequence, the ``LinkerPrimerSequence`` in the mapping file is not important. In the following example, we've simply included an identifier for the LinkerPrimerSequence. An example mapping file to match this data might look like::

	#SampleID	BarcodeSequence	LinkerPrimerSequence	SampleType	Description
	Sample1	AGCTGT	ILBC_17	Freshwater	creek_3-4_cm_depth
	S2	AACTGT	ILBC_19	Ocean	ocean_2-3_cm_depth
	S3	TCATGT	ILBC_21	Ocean	ocean_3-4_cm_depth 
	sample4	GCTGGT	ILBC_26	Feces	fecal_subject_5

Parsing 5' read data
---------------------
To process a single lane of 5' read Illumina data, you would run the command::

	split_libraries_illumina.py -5 s_1_1_sequence.txt -o 5prime/ -m illumina_mapping.txt

Similarly, to process multiple lanes of 5' reads of Illumina data, you would run the command::

	split_libraries_illumina.py -5 s_1_1_sequence.txt,s_2_1_sequence.txt -o 5prime/ -m illumina_mapping.txt
	
Note that the sequence filepaths are comma-separated. There can be no spaces between the filepaths and the comma(s). This would result in four output files in ``5prime/``:
	
	#. ``s_1_5prime_seqs.fasta``: lane 1 sequences which have passed quality filtering
	#. ``s_1_5prime_qual.txt``: lane 1 quality scores for sequences which have passed quality filtering
	#. ``s_2_5prime_seqs.fasta``: lane 2 sequences which have passed quality filtering
	#. ``s_2_5prime_qual.txt``: lane 2 quality scores for sequences which have passed quality filtering

If you wish to combine the data from different lanes for downstream analysis, you can do this using ``cat``. For example::
	
	cat s_1_5prime_seqs.fasta >> s_5prime_seqs.fasta
	cat s_2_5prime_seqs.fasta >> s_5prime_seqs.fasta
	
Parsing 3' read data
---------------------
Parsing 3' read data is handled in the same way as parsing 5' read data, except that the ``-3`` parameter is passed in place of the ``-5`` parameter. For example, to process multiple lanes of 3' Illumina read data, you would run the command::

	split_libraries_illumina.py -3 s_1_2_sequence.txt,s_2_2_sequence.txt -o 3prime/ -m illumina_mapping.txt
	
It is very important that you get the orientation correct, as the sequences are adjusted differently according to whether they are passed as ``-5`` or ``-3``.
	
Parsing paired-end read data
------------------------------
To process a single lane of paired-end read Illumina data, you would run the command::

	split_libraries_illumina.py -5 s_1_1_sequence.txt -3 s_1_2_sequence.txt -o paired_end/ -m illumina_mapping.txt

To process multiple lanes of paired-end read Illumina data, you would run the command::

	split_libraries_illumina.py -5 s_1_1_sequence.txt,s_2_1_sequence.txt -3 s_1_2_sequence.txt,s_2_2_sequence.txt -o paired_end/ -m illumina_mapping.txt

Note that you must correctly match the order of the sequence file paths passed via ``-5`` and ``-3``. In this case, the output sequence files will contain the 5' ends of the reads concatenated with the 3' ends of the reads. The orientations of the reads are adjusted so the full read makes sense. It is important to note however that aside from adjusting the read orientation, nothing is done to modify the read. So, for example, if your reads overlap there will be a repeat at the junction of the reads corresponding to the overlap.

Quality filtering of Illumina data
------------------------------------
The methods currently included in QIIME for quality filtering Illumina data are presented in a publication that has been provisionally accepted. A link to that article will be provided here upon full acceptance. Briefly, a sequence is discarded if any of the following conditions are met:
	
	* The barcode contains one or more ``N`` bases, corresponding to ambiguous base calls.
	* The barcode is not an exact match to a barcode in the mapping file (to disable this, pass ``-u``, which will cause the resulting sequences to be store with sample ID ``Unassigned``.)
	* The sequence contains one or more ``N`` bases, corresponding to ambiguous base calls.
	* The high-quality region of the sequence is less than 75 bases long (adjustable with the ``-p`` parameter), where high-quality regions is defined a stretch of bases containing no more than 1 (adjustable with the ``-r`` parameter) quality score less than 1e-5 (adjustable with the ``-q`` parameter). In other words, with the default parameter settings, the read is truncated at the base preceding the first low quality stretch, and the truncated sequence must be greater than or equal to 75 bases long to be retained. 
	
Using Illumina data in downstream analysis
-------------------------------------------
Once you have processed your Illumina data with `split_libraries_illumina.py <../scripts/split_libraries_illumina.html>`_, you can use the resulting files in all downstream QIIME scripts, including the workflow scripts. Currently the quality files (e.g., ``s_2_5prime_qual.txt``) are not used anywhere in QIIME, but are provided to support more complex quality filtering. We are very interested in collaborating on developing better ways to quality filter Illumina data, so please feel free to get in touch if you have interest in this problem.
