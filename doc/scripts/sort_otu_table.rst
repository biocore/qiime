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
		Input OTU table filepath in BIOM format.
	-o, `-`-output_fp
		Output OTU table filepath.
	
	**[OPTIONAL]**
		
	-m, `-`-mapping_fp
		Input metadata mapping filepath. [default: None]
	-s, `-`-sort_field
		Category to sort OTU table by. [default: None]
	-l, `-`-sorted_sample_ids_fp
		Sorted sample id filepath [default: None]


**Output:**




sort samples by the age field in the mapping file

::

	sort_otu_table.py -i otu_table.biom -o dob_sorted_otu_table.biom -m Fasting_Map.txt -s DOB

sort samples based on order in a file where each line starts with a sample id

::

	sort_otu_table.py -i otu_table.biom -o sorted_otu_table.biom -l sample_id_list.txt


