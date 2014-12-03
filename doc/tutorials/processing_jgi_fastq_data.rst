.. _processing_jgi_fastq_data:

===============================================================================================
Processing Joint Genome Institute fastq files for compatibility with `split_libraries_fastq.py`
===============================================================================================

Introduction
------------

The Joint Genome Institute (JGI) delivers a single fastq file to their users, which is not compatible with `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_. To fix this issue, first we have to split this file into forward and reverse reads and then extract the barcode from the forward file.

Processing a single fastq with reverse and reverse reads
--------------------------------------------------------

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

To split the input file in forward and revers we will use `extract_reads_from_interleaved_file.py <../scripts/extract_reads_from_interleaved_file.html>`_.

	extract_reads_from_interleaved_file.py -i reads_to_extract.fastq -o extracted_reads

Then the next step is to create the barcode (index) fastq file by using the forward read generated in the previous step and `extract_barcodes.py <../scripts/extract_barcodes.html>`_. Note how we are using the parameter `-l 11` because the primers in this example have that length, see header in the example:

	extract_barcodes.py -f extracted_reads/forward_reads.fastq -o barcode_reads -c barcode_in_label --char_delineator '1:N:0:' -l 11
