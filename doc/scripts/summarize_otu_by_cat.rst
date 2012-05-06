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
		Summarize OTU table using this category. The user can also combine columns in the mapping file by separating the categories by "&&" without spaces. [REQUIRED]
	-o, `-`-output_fp
		Output OTU table filepath. [REQUIRED]
	
	**[OPTIONAL]**
		
	-n, `-`-normalize
		Normalize OTU counts, where the OTU table columns sum to 1.


**Output:**




**Example:**

 Collapsed otu_table.biom on the 'Treatment' column in Fasting_Map.txt and write the resulting OTU table to otu_table_by_treatment.txt

::

	summarize_otu_by_cat.py -c otu_table.biom -i Fasting_Map.txt -m Treatment -o otu_table_by_treatment.biom

 Combine two categories and collapse otu_table.biom on the 'Sex' and 'Age' columns in map.txt and write the resulting OTU table to otu_table_by_sex_and_age.txt

::

	summarize_otu_by_cat.py -c otu_table.biom -i Fasting_Map.txt -m 'Treatment&&DOB' -o otu_table_by_treatment_and_dob.biom


