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

The result of `multiple_rarefactions.py <./multiple_rarefactions.html>`_ consists of a number of files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.txt, and are named in the following way: rarefaction_100_0.txt, where "100" corresponds to the sequences per sample and "0" the iteration.


**Examples:**

An example of this script, where the user sets the minimum ("-m") and maximum ("-x") number of sequences per sample to 100 and 1200, respectively, while using steps ("-s") of 100, performing 2 iterations at each sampling depth ("-n"), and outputting the results to the directory "rarefaction_tables/" is shown by the following command:

::

	multiple_rarefactions.py -i otu_table.txt -m 100 -x 1200 -s 100 -n 2 -o rarefaction_tables/

As a result, this command produces subsamples of the input otu_table.txt at 100 seqs per sample (twice), 200 seqs per sample (twice) ... 1200 seqs per sample (twice), which produces 24 rarefied otu talbes in the "rarefaction_tables" directory.

Any sample containing fewer sequences in the input file than the requested number of sequences per sample is removed from the output rarefied otu table. To include samples with fewer than the requested number, you must manually add those samples to the resulting otu tables


