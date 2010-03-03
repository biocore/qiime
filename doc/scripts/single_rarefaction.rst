.. _single_rarefaction:

.. index:: single_rarefaction.py

*single_rarefaction.py* -- Perform rarefaction on a single otu table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, an OTU table.  This does not provide curves of diversity by number of sequences in a sample. Rather it creates a subsampled OTU table by random sampling (without replacement) of the input OTU table.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG.


**Usage:** :file:`single_rarefaction.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input otu table filepath
	-o, `-`-output_path
		Write output rarefied otu tables to this filepath
	-d, `-`-depth
		Sequences per sample to subsample
	
	**[OPTIONAL]**
		
	`-`-small_included
		Samples containing fewer seqs than the rarefaction level are included in the output but not rarefied [default: False]
	`-`-lineages_included
		Output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: False]


**Output:**

The results of `single_rarefaction.py <./single_rarefaction.html>`_ consist of a single subsampled OTU table. The file has the same otu table format as the input otu_table.txt. note: if the output file would be empty, no file is written


**Example:**

subsample otu_table.txt at 400 seqs/sample (-d), write results to a file (i.e. rarefaction_400_17.txt) (samples which have fewer than 400 sequences will be included without subsampleing, exactly as they appear in otu_table.txt

::

	single_rarefaction.py -i otu_table.txt -o rarefaction_400_17.txt -d 400 --small_included

(naming convention implies that the depth is 200 seqs/sam, iteration 17 at that depth (18th file written, due to iter 0))


