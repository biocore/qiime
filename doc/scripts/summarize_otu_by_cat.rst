.. _summarize_otu_by_cat:

.. index:: summarize_otu_by_cat.py

*summarize_otu_by_cat.py* -- Create a summarized OTU table for a specific metadata category
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script generates an otu table where the SampleIDs are replaced by a specific category from the user-generated mapping file. The script uses the OTU file otus.txt (-c) and the user mapping file meta.txt. The user must also specify a metadata category equivalent to one of the column names in the mapping file. If the user wants the counts to be normalized by sample use th normalize flag (-n) the default is False meaning it is only the raw counts. The output is a file called <meta category>_otu_table.txt, it will be put int the current working directory unless specified by the user (-o).


**Usage:** :file:`summarize_otu_by_cat.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_map
		Name of input map file [REQUIRED]
	-c, `-`-otu_file
		Name of otu table file [REQUIRED]
	-m, `-`-meta_category
		Name of category for OTU table [REQUIRED]
	
	**[OPTIONAL]**
		
	-o, `-`-dir-prefix
		Directory prefix for all analyses [default: cwd]
	-n, `-`-normalize_flag
		If True will normalize counts [default: False]


**Output:**

The output is an otu table called <meta category>_otu_table.txt, 


**Example:**

Create an otu table for a user specified category. This script uses an OTU table (otu_table.txt) and a user-generated mapping file (mapping_file.txt). The user must also specify a metadata category equivalent to one of the column names in their mapping file (i.e. time). If the user wants the counts to be normalized by sample, they can use the normalize flag (-n), however; the default value for this flag is False, which means it will use the raw counts. The resulting files will be it will be written in the current working directory, unless specified by the user (-o).

::

	summarize_otu_by_cat.py -c otu_table.txt -i mapping_file.txt -m time -o qiime_run/ -n


