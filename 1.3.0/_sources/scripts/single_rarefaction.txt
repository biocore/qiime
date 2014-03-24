.. _single_rarefaction:

.. index:: single_rarefaction.py

*single_rarefaction.py* -- Perform rarefaction on an otu table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, an OTU table.  This does not provide curves of diversity by number of sequences in a sample. Rather it creates a subsampled OTU table by random sampling (without replacement) of the input OTU table.  Samples that have fewer sequences then the requested rarefaction depth are omitted from the ouput otu tables.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG.


**Usage:** :file:`single_rarefaction.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input OTU table filepath.
	-o, `-`-output_path
		Output OTU table filepath.
	-d, `-`-depth
		Number of sequences to subsample per sample.
	
	**[OPTIONAL]**
		
	`-`-lineages_included
		Deprecated: lineages are now included by default. Pass --supress_lineages_included to prevent output OTU tables from including taxonomic (lineage) information for each OTU. Note: this will only work if lineage information is in the input OTU table.
	`-`-suppress_lineages_included
		Exclude taxonomic (lineage) information for each OTU.
	-k, `-`-keep_empty_otus
		Retain OTUs of all zeros, which are usually omitted from the output OTU tables. [default: False]


**Output:**

The results of `single_rarefaction.py <./single_rarefaction.html>`_ consist of a single subsampled OTU table. The file has the same otu table format as the input otu_table.txt. note: if the output file would be empty, no file is written


**Example:**

subsample otu_table.txt at 400 seqs/sample (-d), write results to a file (i.e. rarefaction_400_17.txt) 

::

	single_rarefaction.py -i otu_table.txt -o rarefaction_400_17.txt -d 400

(naming convention rarefaction_400_17.txt implies that the depth is 400 seqs/sam, iteration 17 at that depth (18th file written, due to iter 0))


