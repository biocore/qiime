.. _compute_core_microbiome:

.. index:: compute_core_microbiome.py

*compute_core_microbiome.py* -- Identify the core microbiome.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`compute_core_microbiome.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input otu table in BIOM format
	-o, `-`-output_dir
		Directory to store output data
	
	**[OPTIONAL]**
		
	`-`-max_fraction_for_core
		The max fractions of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: 1.0]
	`-`-min_fraction_for_core
		The min fractions of samples that an OTU must be observed in to be considered part of the core as a number in the range [0,1] [default: 0.5]
	`-`-num_fraction_for_core_steps
		The number of evenly sizes steps to take between min_fraction_for_core and max_fraction_for_core [default: 11]
	`-`-otu_md
		The otu metadata category to write to the output file [defualt: taxonomy]
	`-`-mapping_fp
		Mapping file path (for use with --valid_states) [default: None]
	`-`-valid_states
		Description of sample ids to retain (for use with --mapping_fp) [default: None]


**Output:**




Identify the core OTUs in otu_table.biom, defined as the OTUs that are present in at least 50% of the samples. Write the list of core OTUs to a text file, and a new BIOM file containing only the core OTUs.

::

	compute_core_microbiome.py -i otu_table.biom -o otu_table_core

Identify the core OTUs in otu_table.biom, defined as the OTUs that are present in all of the samples in the 'Fast' treatment (as specified in the mapping file). Write the list of core OTUs to a text file.

::

	compute_core_microbiome.py -i otu_table.biom -o otu_table_core_fast --mapping_fp map.txt --valid_states "Treatment:Fast"


