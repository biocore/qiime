.. _process_sff:

.. index:: process_sff.py

*process_sff.py* -- Convert sff to FASTA and QUAL files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script converts a directory of sff files into FASTA and QUAL files.

This script requires that 454's off-instrument apps (sffinfo, sfffile) are in your path.


**Usage:** :file:`process_sff.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Input directory of sff files


**Output:**

This script results in FASTA and QUAL formatted files.


**Simple example:**

Convert all the sffs in directory "sffs/" to fasta and qual.

::

	process_sff.py -i sffs/


