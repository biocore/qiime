.. _multiple_rarefactions:

.. index:: multiple_rarefactions.py

*multiple_rarefactions.py* -- Perform multiple subsamplings/rarefactions on an otu table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To perform bootstrap, jackknife, and rarefaction analyses, the otu table must be subsampled (rarefied).  This script rarefies, or subsamples, OTU tables.  This does not provide curves of diversity by number of sequences in a sample.  Rather it creates a series of subsampled OTU tables by random sampling (without replacement) of the input OTU table.  Samples that have fewer sequences then the requested rarefaction depth for a given output otu table are omitted from those ouput otu tables.  The pseudo-random number generator used for rarefaction by subsampling is NumPy's default - an implementation of the Mersenne twister PRNG.


**Usage:** :file:`multiple_rarefactions.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input OTU table filepath.
	-o, `-`-output_path
		Output directory.
	-m, `-`-min
		Minimum number of seqs/sample for rarefaction.
	-x, `-`-max
		Maximum number of seqs/sample (inclusive) for rarefaction. 
	-s, `-`-step
		Size of each steps between the min/max of seqs/sample (e.g. min, min+step... for level <= max).
	
	**[OPTIONAL]**
		
	-n, `-`-num-reps
		The number of iterations at each step. [default: 10]
	`-`-lineages_included
		Retain taxonomic (lineage) information for each OTU. Note: this will only work if lineage information is in the input OTU table. [default: False]
	-k, `-`-keep_empty_otus
		Retain OTUs of all zeros, which are usually omitted from the output OTU tables. [default: False]


**Output:**

The result of `multiple_rarefactions.py <./multiple_rarefactions.html>`_ consists of a number of biom files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.biom, and are named in the following way: rarefaction_100_0.biom, where "100" corresponds to the sequences per sample and "0" the iteration.


**Generate rarefied OTU tables:**

Generate rarefied OTU tables beginning with 10 (-m) sequences/sample through 140 (-x) sequences per sample in steps of of 10 (-s), performing 2 iterations at each sampling depth (-n). All resulting OTU tables will be written to 'rarefied_otu_tables' (-o). Any sample containing fewer sequences in the input file than the requested number of sequences per sample is removed from the output rarefied otu table.

::

	multiple_rarefactions.py -i otu_table.biom -m 10 -x 140 -s 10 -n 2 -o rarefied_otu_tables/


