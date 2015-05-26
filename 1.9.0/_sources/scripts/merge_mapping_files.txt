.. _merge_mapping_files:

.. index:: merge_mapping_files.py

*merge_mapping_files.py* -- Merge mapping files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script provides a convenient interface for merging mapping files which contain data on different samples.


**Usage:** :file:`merge_mapping_files.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fps
		The input mapping files in a comma-separated list
	-o, `-`-output_fp
		The output mapping file to write
	
	**[OPTIONAL]**
		
	-n, `-`-no_data_value
		Value to represent missing data (i.e., when all fields are not defined in all mapping files) [default: no_data]
	`-`-case_insensitive
		If present the headers will be merged case insensitivly and transformed to upper case [default: False]


**Output:**

The result of this script is a merged mapping file (tab-delimited).


**Example:**

Merge two mapping files into a new mapping file (merged_mapping.txt). In cases where a mapping field is not provided for some samples, add the value 'Data not collected'.

::

	merge_mapping_files.py -m map_controls.txt,map_fasting.txt -o merged_mapping.txt -n 'Data not collected'


