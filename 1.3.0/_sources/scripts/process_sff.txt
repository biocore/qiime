.. _process_sff:

.. index:: process_sff.py

*process_sff.py* -- Convert sff to FASTA and QUAL files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script converts a directory of sff files into FASTA, QUAL and flowgram files.



**Usage:** :file:`process_sff.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Input directory of sff files
	
	**[OPTIONAL]**
		
	-f, `-`-make_flowgram
		Generate a flowgram file. [default: False]
	-t, `-`-convert_to_FLX
		Convert Titanium reads to FLX length. [default: False]
	`-`-use_sfftools
		Use the external programs sfffile and sffinfo for processing, instead of the equivalent python implementation
	-o, `-`-output_dir
		Input directory of sff files [default: same as input dir]


**Output:**

This script results in FASTA and QUAL formatted files.


**Simple example:**

Convert all the sffs in directory "sffs/" to fasta and qual.

::

	process_sff.py -i sffs/

**Flowgram example:**

Convert all the sffs in directory "sffs/" to fasta and qual, along with a flowgram file.

::

	process_sff.py -i sffs/ -f

**Output example:**

Convert all the sffs in directory "sffs/" to fasta and qual, along with a flowgram file and write them to another directory.

::

	process_sff.py -i sffs/ -f -o output_dir


