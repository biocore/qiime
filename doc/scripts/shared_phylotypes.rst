.. _shared_phylotypes:

.. index:: shared_phylotypes.py

*shared_phylotypes.py* -- Compute shared OTUs between all pairs of samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script computes from an OTU table a matrix with the number of shared phylotypes between all pairs of samples.


**Usage:** :file:`shared_phylotypes.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table in biom format or a directory containing OTU tables
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-r, `-`-reference_sample
		Name of reference sample to which all pairs of samples should be compared [default: None]


**Output:**




**Single example:**

Compute shared OTUs on one OTU table for all samples

::

	shared_phylotypes.py -i otu_table.biom -o shared_otus.txt

**Reference sample example:**

Compute shared OTUs with respect to a reference sample. Computes shared OTUs between all pairs of samples and the reference sample. E.g. in a transplant study this can be used to establish a base line count of shared OTUs with the Donor sample before and after the transplant.

::

	shared_phylotypes.py -i otu_table.biom -o shared_otus_PC.636.txt -r PC.636

**Batch mode example:**

Compute shared OTUs for a set of OTU tables, e.g. from running `multiple_rarefactions.py <./multiple_rarefactions.html>`_, with an even number of sequences per sample. The resulting directory can be fed to `dissimilarity_mtx_stats.py <./dissimilarity_mtx_stats.html>`_, which computes mean, median and the standard deviation on the provided tables.

::

	shared_phylotypes.py -i rarefied_otu_tables/ -o shared_otus/


