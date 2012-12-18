.. _filter_otus_from_otu_table:

.. index:: filter_otus_from_otu_table.py

*filter_otus_from_otu_table.py* -- Filter OTUs from an OTU table based on their observation counts or identifier.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`filter_otus_from_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input otu table filepath in biom format
	-o, `-`-output_fp
		The output filepath in biom format
	
	**[OPTIONAL]**
		
	`-`-negate_ids_to_exclude
		Keep OTUs in otu_ids_to_exclude_fp rather than discard them [default:False] 
	-n, `-`-min_count
		The minimum total observation count of an otu for that otu to be retained [default: 0]
	`-`-min_count_fraction
		Fraction of the total observation (sequence) count to apply as the minimum total observation count of an otu for that otu to be retained. this is a fraction, not percent, so if you want to filter to 1%, you specify 0.01. [default: 0]
	-x, `-`-max_count
		The maximum total observation count of an otu for that otu to be retained [default: infinity]
	-s, `-`-min_samples
		The minimum number of samples an OTU must be observed in for that otu to be retained [default: 0]
	-y, `-`-max_samples
		The maximum number of samples an OTU must be observed in for that otu to be retained [default: infinity]
	-e, `-`-otu_ids_to_exclude_fp
		File containing list of OTU ids to exclude: can be a text file with one id per line, a text file where id is the first value in a tab-separated line, or can be a fasta file (extension must be .fna or .fasta) [default: None]


**Output:**




**Singleton filtering:**

Discard all OTUs that are observed fewer than 2 times (i.e., singletons)

::

	filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_no_singletons.biom -n 2

**Abundance filtering:**

Discard all OTUs that are observed greater than 100 times (e.g., if you want to look at low abundance OTUs only)

::

	filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_low_abundance.biom -x 100

**Chimera filtering:**

Discard all OTUs listed in chimeric_otus.txt (e.g., to remove chimeric OTUs from an OTU table)

::

	filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_non_chimeric.biom -e chimeric_otus.txt


