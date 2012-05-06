.. _trflp_file_to_otu_table:

.. index:: trflp_file_to_otu_table.py

*trflp_file_to_otu_table.py* -- Convert TRFLP text file to an OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The input for this script is a TRLFP text file. The output of this script is an OTU table text file that can be use with QIIME for further analysis 


**Usage:** :file:`trflp_file_to_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path: TRFLP text file
	-o, `-`-output_path
		Output file: OTU table


**Output:**




**Usage:**

You need to pass a TRFLP text, the script will remove not wanted chars sample and otus names, and will add zeros as need it

::

	trflp_file_to_otu_table.py -i trflp_in.txt -o otu_table.biom


