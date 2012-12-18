.. _filter_taxa_from_otu_table:

.. index:: filter_taxa_from_otu_table.py

*filter_taxa_from_otu_table.py* -- Filter taxa from an OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This scripts filters an OTU table based on taxonomic metadata. It can be applied for positive filtering (i.e., keeping only certain taxa), negative filtering (i.e., discarding only certain taxa), or both at the same time.


**Usage:** :file:`filter_taxa_from_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_otu_table_fp
		The input otu table filepath
	-o, `-`-output_otu_table_fp
		The output otu table filepath
	
	**[OPTIONAL]**
		
	-p, `-`-positive_taxa
		Comma-separated list of taxa to retain [default: None; retain all taxa]
	-n, `-`-negative_taxa
		Comma-separated list of taxa to discard [default: None; retain all taxa]
	`-`-metadata_field
		Observation metadata identifier to filter based on [default: taxonomy]


**Output:**




Filter otu_table.biom to include only OTUs identified as __Bacteroidetes or p__Firmicutes.

::

	filter_taxa_from_otu_table.py -i otu_table.biom -o otu_table_bac_firm_only.biom -p p__Bacteroidetes,p__Firmicutes

Filter otu_table.biom to exclude OTUs identified as p__Bacteroidetes or p__Firmicutes.

::

	filter_taxa_from_otu_table.py -i otu_table.biom -o otu_table_non_bac_firm.biom -n p__Bacteroidetes,p__Firmicutes

Filter otu_table.biom to include OTUs identified as p__Firmicutes but not c__Clostridia.

::

	filter_taxa_from_otu_table.py -i otu_table.biom -o otu_table_all_firm_but_not_clos.biom -p p__Firmicutes -n c__Clostridia


