.. _process_qseq:

.. index:: process_qseq.py

*process_qseq.py* -- Given a directory of per-swath qseq files, this script generates a single fastq per lane.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`process_qseq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		The input directory
	-o, `-`-output_dir
		The output directory
	-r, `-`-read
		The read number to consider
	
	**[OPTIONAL]**
		
	-l, `-`-lanes
		The lane numbers to consider, comma-separated [defaut: 1,2,3,4,5,6,7,8]
	-b, `-`-bases
		The number of bases to include (useful for slicing a barcode) [defaut: all]
	`-`-ignore_pass_filter
		Ignore the illumina pass filter [default:False; reads with 0 in pass  filter field are discarded]


**Output:**




Generate fastq files from all lanes of read 1 data in the current directory.

::

	process_qseq.py -i ./ -o ./fastq/ -r 1

Generate fastq files from all lanes of read 2 data in the current directory, truncating the sequences after the first 12 bases.

::

	process_qseq.py -i ./ -o ./fastq/ -r 2 -b 12


