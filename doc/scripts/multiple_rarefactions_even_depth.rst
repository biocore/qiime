.. _multiple_rarefactions_even_depth:

.. index:: multiple_rarefactions_even_depth.py

*multiple_rarefactions_even_depth.py* -- Perform multiple rarefactions on a single otu table, at one depth of sequences/sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, an OTU table.  This does not provide curves of diversity by number of sequences in a sample.  Rather it creates a subsampled OTU table by random sampling (without replacement) of the input OTU table.  Samples that have fewer sequences then the requested rarefaction depth are omitted from the ouput otu tables.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG.


**Usage:** :file:`multiple_rarefactions_even_depth.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input otu table filepath
	-o, `-`-output_path
		Write output rarefied otu tables files to this dir (makes dir if it doesn't exist)
	-d, `-`-depth
		Sequences per sample to subsample
	
	**[OPTIONAL]**
		
	-n, `-`-num-reps
		Num iterations at each seqs/sample level [default: 10]
	`-`-lineages_included
		Output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: False]
	-k, `-`-keep_empty_otus
		Otus (rows) of all zeros are usually omitted from the output otu tables, with -k they will not be removed from the output files [default: False]


**Output:**

The results of this script consist of n subsampled OTU tables, written to the directory specified by -o. The file has the same otu table format as the input otu_table.biom. Note: if the output files would be empty, no files are written.


**Example:**

subsample otu_table.biom at 100 seqs/sample (-d) 10 times (-n) and write results to files (e.g., rarefaction_400_0.biom) in 'rarefied_otu_tables/' (-o).

::

	multiple_rarefactions_even_depth.py -i otu_table.biom -o rarefied_otu_tables/ -d 100 -n 10


