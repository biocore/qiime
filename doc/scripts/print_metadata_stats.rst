.. _print_metadata_stats:

.. index:: print_metadata_stats.py

*print_metadata_stats.py* -- Count the number of samples associated to a category value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Sum up the number of samples with each category value and print this information.


**Usage:** :file:`print_metadata_stats.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_file
		The input metadata file
	-c, `-`-category
		The category to examine
	
	**[OPTIONAL]**
		
	-o, `-`-output_fp
		Path where output will be written [default: print to screen]


**Output:**

Two columns, the first being the category value and the second being the count. Output is to standard out. If there are unspecified values, the output category is identified as ***UNSPECIFIED***


**Example:**

Count the number of samples associated with Treatment

::

	print_metadata_stats.py -m $PWD/mapping.txt -c Treatment

**Example writting the output to a file:**

Count the number of samples associated with Treatment and save them to a file called stats.txt

::

	print_metadata_stats.py -m mapping.txt -c Treatment -o stats.txt


