.. _pick_otus_through_otu_table:

.. index:: pick_otus_through_otu_table.py

*pick_otus_through_otu_table.py* -- A workflow script for picking OTUs through building OTU tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a sequence file and performs all processing steps through building the OTU table.


**Usage:** :file:`pick_otus_through_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input fasta file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]


**Output:**

This script will produce an OTU mapping file (`pick_otus.py <./pick_otus.html>`_), a representative set of sequences (FASTA file from `pick_rep_set.py <./pick_rep_set.html>`_), a sequence alignment file (FASTA file from `align_seqs.py <./align_seqs.html>`_), taxonomy assignment file (from `assign_taxonomy.py <./assign_taxonomy.html>`_), a filtered sequence alignment (from `filter_alignment.py <./filter_alignment.html>`_), a phylogenetic tree (Newick file from `make_phylogeny.py <./make_phylogeny.html>`_) and an OTU table (from `make_otu_table.py <./make_otu_table.html>`_).


**Simple example:**

The following command will start an analysis on inseq1.fasta (-i), which is a post-split_libraries fasta file. The sequence identifiers in this file should be of the form <sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data preparation, are applied.

1. Pick OTUs with uclust at similarity of 0.97;

2. Pick a representative set with the most_abundant method;

3. Align the representative set with PyNAST;

4. Assign taxonomy with RDP classifier;

5. Filter the alignment prior to tree building - remove positions which are all gaps, and specified as 0 in the lanemask;

6. Build a phylogenetic tree with FastTree;

7. Build an OTU table.

All output files will be written to the directory specified by -o, and 
subdirectories as appropriate.


::

	pick_otus_through_otu_table.py -i inseqs1.fasta -o wf1/ -p custom_parameters.txt


