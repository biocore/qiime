.. _summarize_otu_by_cat:

.. index:: summarize_otu_by_cat.py

*summarize_otu_by_cat.py* -- Summarize an OTU table by a single column in the mapping file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Collapse an OTU table based on values in a single column in the mapping file. For example, if you have 10 samples, five of which are from females and five of which are from males, you could use this script to collapse the ten samples into two corresponding based on their values in a 'Sex' column in your mapping file.


**Usage:** :file:`summarize_otu_by_cat.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-mapping_fp
		Input metadata mapping filepath [REQUIRED]
	-c, `-`-otu_table_fp
		Input OTU table filepath. [REQUIRED]
	-m, `-`-mapping_category
		Summarize OTU table using this category. [REQUIRED]
	-o, `-`-output_fp
		Output OTU table filepath. [REQUIRED]
	
	**[OPTIONAL]**
		
	-n, `-`-normalize
		Normalize OTU counts, where the OTU table columns sum to 1.


**Output:**




**Example:**

 Collapsed otu_table.txt on the 'Sex' column in map.txt and write the resulting OTU table to otu_table_by_sex.txt

::

	summarize_otu_by_cat.py -c otu_table.txt -i map.txt -m Sex -o otu_table_by_sex.txt


