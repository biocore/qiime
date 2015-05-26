.. _insert_seqs_into_tree:

.. index:: insert_seqs_into_tree.py

*insert_seqs_into_tree.py* -- Tree Insertion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a set of aligned sequences (query) either in the same file as the aligned reference set or separated (depending on method) along with a starting tree and produces a new tree containing the query sequences. This script requires that the user is running Raxml v7.3.0, PPlacer git repository version and ParsInsert 1.0.4.


**Usage:** :file:`insert_seqs_into_tree.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-o, `-`-output_dir
		Path to the output directory
	-t, `-`-starting_tree_fp
		Starting Tree which you would like to insert into.
	-r, `-`-refseq_fp
		Filepath for reference alignment
	
	**[OPTIONAL]**
		
	-m, `-`-insertion_method
		Method for aligning sequences. Valid choices are: pplacer, raxml_v730, parsinsert [default: raxml_v730]
	-s, `-`-stats_fp
		Stats file produced by tree-building software. REQUIRED if -m pplacer [default: None]
	-p, `-`-method_params_fp
		Parameters file containing method-specific parameters to use. [default: None]


**Output:**

The result of this script produces a tree file (in Newick format) along with a log file containing the output from the underlying tool used for tree insertion.


**RAxML Example (default):**

If you just want to use the default options, you can supply an alignment files where the query and reference sequences are included, along with a starting tree as follows:

::

	insert_seqs_into_tree.py -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results

**ParsInsert Example:**

If you want to insert sequences using pplacer, you can supply a fasta file containg query sequences (aligned to reference sequences) along with the reference alignment, a starting tree and the stats file produced when building the starting tree via pplacer as follows:

::

	insert_seqs_into_tree.py -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results -m parsinsert

**Pplacer Example:**

If you want to insert sequences using pplacer, you can supply a fasta file containg query sequences (aligned to reference sequences) along with the reference alignment, a starting tree and the stats file produced when building the starting tree via pplacer as follows:

::

	insert_seqs_into_tree.py -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results -m pplacer

**Parameters file:**

Additionally, users can supply a parameters file to change the options of the underlying tools as follows:

::

	insert_seqs_into_tree.py -i aligned_query_seqs.fasta -r aligned_reference_seqs.fasta -t starting_tree.tre -o insertion_results -p raxml_parameters.txt


