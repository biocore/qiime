.. _split_otu_table:

.. index:: split_otu_table.py

*split_otu_table.py* -- Split in a single OTU table into one OTU table per value in a specified field of the mapping file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`split_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input otu table
	-m, `-`-mapping_fp
		The mapping file path
	-f, `-`-mapping_field
		Mapping column to split otu table on
	-o, `-`-output_dir
		The output directory


**Output:**




Split otu_table.biom into per-study OTU tables, and store the results in ./per_study_otu_tables/

::

	split_otu_table.py -i otu_table.biom -m Fasting_Map.txt -f Treatment -o per_study_otu_tables


