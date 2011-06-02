.. _per_library_stats:

.. index:: per_library_stats.py

*per_library_stats.py* -- Calculate per library statistics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Given an otu table, compute and print the (min, max, median, mean) number of seqs per library.


**Usage:** :file:`per_library_stats.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	
	**[OPTIONAL]**
		
	-m, `-`-mapfile
		A mapping file. If included, this script will modify the mapping file to include sequences per sample (library) information, and write the modified mapping file to the path specified by -o. The sequences (individuals) per sample is presented in a new column entitled "NumIndividuals", and samples present in the mapping file but not the otu table have the value "na" in this column. Note also that the location of comments is not preserved in the new mapping file
	-o, `-`-outputfile
		The output filepath where the modified mapping file will be written


**Output:**

The resulting statistics are written to stdout. If -m is passed, a new mapping file is written to the path specified by -o, in addition to the statistics written to stdout


**Example:**

Calculate statistics on an OTU table (otu_table.txt)

::

	per_library_stats.py -i otu_table.txt

**Example appending results to mapping file:**

Calculate statistics on an OTU table (otu_table.txt)

::

	per_library_stats.py -i otu_table.txt -m old_map.txt -o new_map.txt


