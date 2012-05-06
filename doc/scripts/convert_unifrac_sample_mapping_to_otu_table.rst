.. _convert_unifrac_sample_mapping_to_otu_table:

.. index:: convert_unifrac_sample_mapping_to_otu_table.py

*convert_unifrac_sample_mapping_to_otu_table.py* -- Convert a UniFrac sample mapping file to an OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script allows users that have already created sample mapping (environment) files for use with the Unifrac web interface to use QIIME. QIIME records this data in an OTU table.


**Usage:** :file:`convert_unifrac_sample_mapping_to_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-sample_mapping_fp
		Path to the sample mapping file
	-o, `-`-output_fp
		Path to output file


**Output:**

The result of this script is an OTU table.


**Example:**

Convert a UniFrac sample mapping (environment) file into a biom-formatted OTU table: 

::

	convert_unifrac_sample_mapping_to_otu_table.py -i otu_table.sample_mapping.txt -o otu_table.biom


