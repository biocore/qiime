.. _split_otu_table_by_taxonomy:

.. index:: split_otu_table_by_taxonomy.py

*split_otu_table_by_taxonomy.py* -- Script to split a single OTU table into multiple tables based on the taxonomy at some user-specified depth.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`split_otu_table_by_taxonomy.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input otu table
	-L, `-`-level
		The level to split on
	-o, `-`-output_dir
		The output directory


**Output:**




Split seqs_otu_table.txt into taxon-specific OTU tables based on the third level in the taxonomy, and write the taxon-specific OTU tables to ./L3/

::

	split_otu_table_by_taxonomy.py -i seqs_otu_table.txt -L 3 -o ./L3/


