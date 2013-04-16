.. _pick_de_novo_otus:

.. index:: pick_de_novo_otus.py

*pick_de_novo_otus.py* -- A workflow for de novo OTU picking, taxonomy assignment, phylogenetic tree construction, and OTU table construction.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a sequence file and performs all processing steps through building the OTU table.


**Usage:** :file:`pick_de_novo_otus.py [options]`

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
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 2]


**Output:**

This script will produce an OTU mapping file (`pick_otus.py <./pick_otus.html>`_), a representative set of sequences (FASTA file from `pick_rep_set.py <./pick_rep_set.html>`_), a sequence alignment file (FASTA file from `align_seqs.py <./align_seqs.html>`_), taxonomy assignment file (from `assign_taxonomy.py <./assign_taxonomy.html>`_), a filtered sequence alignment (from `filter_alignment.py <./filter_alignment.html>`_), a phylogenetic tree (Newick file from `make_phylogeny.py <./make_phylogeny.html>`_) and a biom-formatted OTU table (from `make_otu_table.py <./make_otu_table.html>`_).


**Simple example:**

The following command will start an analysis on seqs.fna (-i), which is a post-split_libraries fasta file. The sequence identifiers in this file should be of the form <sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data preparation, are applied: Pick de novo OTUs at 97%; pick a representative sequence for each OTU (the OTU centroid sequence); align the representative set with PyNAST; assign taxonomy with RDP classifier; filter the alignment prior to tree building - remove positions which are all gaps, and specified as 0 in the lanemask; build a phylogenetic tree with FastTree; build an OTU table. All output files will be written to the directory specified by -o, and subdirectories as appropriate. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/otus/


