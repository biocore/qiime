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
	
	**[OPTIONAL]**
		
	-c, `-`-column_rename_ids
		Mapping column used as sample id in the output files. Has to be unique in the splited samples. This option can be helpful to create otu tables and mapping files for Procustes analysis.
	`-`-include_repeat_cols
		By default the new mapping files will not have the columns that have the same information, to include them use this option. This can be helpful to create mapping files for Procrustes analysis.


**Output:**




Split otu_table.txt into per-study OTU tables, and store the results in ./per_study_otu_tables/

::

	split_otu_table.py -i otu_table.txt -m mapping.txt -f Study -o per_study_otu_tables


