.. _convert_otu_table_to_unifrac_sample_mapping:

.. index:: convert_otu_table_to_unifrac_sample_mapping.py

*convert_otu_table_to_unifrac_sample_mapping.py* -- Convert a QIIME OTU table to a UniFrac sample mapping file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script allows users who have picked OTUs in QIIME to convert it to a sample mapping (environment) file for use with the Unifrac web interface.


**Usage:** :file:`convert_otu_table_to_unifrac_sample_mapping.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the otu table
	-o, `-`-output_fp
		Path to output file


**Output:**

The result of this script is a sample mapping file for the UniFrac web interface.


**Example:**

Convert an OTU table (e.g. otu_table.txt) to a unifrac sample mapping (environment) file: 

::

	convert_otu_table_to_unifrac_sample_mapping.py -i otu_table.txt -o sample_mapping.txt


