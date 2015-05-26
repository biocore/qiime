.. _filter_samples_from_otu_table:

.. index:: filter_samples_from_otu_table.py

*filter_samples_from_otu_table.py* -- Filters samples from an OTU table on the basis of the number of observations in that sample, or on the basis of sample metadata. Mapping file can also be filtered to the resulting set of sample ids.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`filter_samples_from_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input otu table filepath in biom format
	-o, `-`-output_fp
		The output filepath in biom format
	
	**[OPTIONAL]**
		
	-m, `-`-mapping_fp
		Path to the map file [default: None]
	`-`-output_mapping_fp
		Path to write filtered mapping file [default: filtered mapping file is not written]
	`-`-sample_id_fp
		Path to file listing sample ids to keep [default: None]
	-s, `-`-valid_states
		String describing valid states (e.g. 'Treatment:Fasting') [default: None]
	-n, `-`-min_count
		The minimum total observation count in a sample for that sample to be retained [default: 0]
	-x, `-`-max_count
		The maximum total observation count in a sample for that sample to be retained [default: infinity]


**Output:**




**Abundance filtering (low coverage):**

Filter samples with fewer than 150 observations from the otu table.

::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_no_low_coverage_samples.biom -n 150

**Abundance filtering (high coverage):**

Filter samples with greater than 149 observations from the otu table.

::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_no_high_coverage_samples.biom -x 149

**Metadata-based filtering (positive):**

Filter samples from the table, keeping samples where the value for 'Treatment' in the mapping file is 'Control'

::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_control_only.biom -m map.txt -s 'Treatment:Control'

**Metadata-based filtering (negative):**

Filter samples from the table, keeping samples where the value for 'Treatment' in the mapping file is not 'Control'

::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_not_control.biom -m map.txt -s 'Treatment:*,!Control'

**List-based filtering:**

Filter samples where the id is listed in samples_to_keep.txt

::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_samples_to_keep.biom --sample_id_fp samples_to_keep.txt


