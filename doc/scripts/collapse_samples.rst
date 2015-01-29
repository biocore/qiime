.. _collapse_samples:

.. index:: collapse_samples.py

*collapse_samples.py* -- Collapse samples in a BIOM table and mapping file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Collapse samples in a BIOM table and mapping file. Values in the BIOM table are collapsed in one of several different ways; see the available options for --collapse_mode. Values in the mapping file are collapsed by grouping the values if they differ for the grouped samples, and by providing the single value if they don't differ for the grouped samples.


**Usage:** :file:`collapse_samples.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-b, `-`-input_biom_fp
		The biom table containing the samples to be collapsed
	-m, `-`-mapping_fp
		The sample metdata mapping file
	`-`-output_biom_fp
		Path where collapsed biom table should be written
	`-`-output_mapping_fp
		Path where collapsed mapping file should be written
	`-`-collapse_fields
		Comma-separated list of fields to collapse on
	
	**[OPTIONAL]**
		
	`-`-collapse_mode
		The mechanism for collapsing counts within groups; valid options are: mean, sum, random, median, first
	`-`-normalize
		Normalize observation counts to relative abundances, so the counts within each sample sum to 1.0. [default: False]


**Output:**

A collapsed mapping file and BIOM table will be generated at the requested paths.


**Collapse samples in biom table and mapping file:**

Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for subject in the mapping file.

::

	collapse_samples.py -b table.biom -m map.txt --output_biom_fp collapsed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields subject

**Collapse samples in biom table and mapping file:**

Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for both subject and replicate-group in the mapping file.

::

	collapse_samples.py -b table.biom -m map.txt --output_biom_fp collapsed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields replicate-group,subject

**Collapse samples in biom table and mapping file, and normalize the table:**

Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for both subject and replicate-group in the mapping file. Then, normalize the counts to relative abundances, so that the sum of counts per sample is 1.0.

::

	collapse_samples.py -b table.biom -m map.txt --output_biom_fp collapsed-normed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields replicate-group,subject --normalize


