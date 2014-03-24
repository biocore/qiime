.. _merge_otu_tables:

.. index:: merge_otu_tables.py

*merge_otu_tables.py* -- Merge two or more OTU tables into a single OTU table.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script merges two or more OTU tables into a single OTU table. This is useful, for example, when you've created several reference-based OTU tables for different analyses and need to combine them for a larger analysis. 

Requirements: It is also very important that your OTUs are consistent across the different OTU tables. For example, you cannot safely merge OTU tables from two independent de novo OTU picking runs. Finally, either all or none of the OTU tables can contain taxonomic information: you can't merge some OTU tables with taxonomic data and some without taxonomic data.


**Usage:** :file:`merge_otu_tables.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		The otu tables in biom format (comma-separated)
	-o, `-`-output_fp
		The output otu table filepath


**Output:**




Merge two OTU tables into a single OTU table

::

	merge_otu_tables.py -i otu_table1.biom,otu_table2.biom -o merged_otu_table.biom


