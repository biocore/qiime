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
		Path to file listing sample ids to keep. Valid formats for the file are: 1) any white space, newline, or tab delimited list of samples, 2) a mapping file with samples in the first column [default: None]
	-s, `-`-valid_states
		String describing valid states (e.g. 'Treatment:Fasting') [default: None]
	-n, `-`-min_count
		The minimum total observation count in a sample for that sample to be retained [default: 0]
	-x, `-`-max_count
		The maximum total observation count in a sample for that sample to be retained [default: infinity]
	`-`-negate_sample_id_fp
		Discard samples specified in --sample_id_fp instead of keeping them [default: False]


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

**ID-based filtering:**

Keep samples where the id is listed in ids.txt

::

	filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom --sample_id_fp ids.txt

**ID-based filtering (negation):**

Discard samples where the id is listed in ids.txt

::

	filter_samples_from_otu_table.py -i otu_table.biom -o filtered_otu_table.biom --sample_id_fp ids.txt --negate_sample_id_fp


