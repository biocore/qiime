.. _sort_otu_table:

.. index:: sort_otu_table.py

*sort_otu_table.py* -- Script for sorting the sample IDs in an OTU table based on a specified value in a mapping file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`sort_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_otu_table
		The input otu table
	-o, `-`-output_fp
		Output otu table filepath
	-m, `-`-mapping_fp
		The mapping file
	-s, `-`-sort_field
		Field to sort by


**Output:**




sort samples by the age field in the mapping file

::

	sort_otu_table.py -i otu_table.txt -o age_sorted_otu_table.txt -m map.txt -s Age


