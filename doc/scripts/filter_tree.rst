.. _filter_tree:

.. index:: filter_tree.py

*filter_tree.py* -- This script prunes a tree based on a set of tip names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a tree and a list of OTU IDs (in one of several supported formats) and outputs a subtree retaining only the tips on the tree which are found in the inputted list of OTUs (or not found, if the --negate option is provided).


**Usage:** :file:`filter_tree.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_tree_filepath
		Input tree filepath
	-o, `-`-output_tree_filepath
		Output tree filepath
	
	**[OPTIONAL]**
		
	-n, `-`-negate
		If negate is True will remove input tips/seqs, if    negate is False, will retain input tips/seqs [default: False]
	-t, `-`-tips_fp
		A list of tips (one tip per line) or sequence identifiers   (tab-delimited lines with a seq identifier in the first field)   which should be retained   [default: None]
	-f, `-`-fasta_fp
		A fasta file where the seq ids should be retained [default: None]


**Output:**

Output is a pruned tree in newick format.


**Prune a tree to include only the tips in tips_to_keep.txt:**

::

	filter_tree.py -i rep_seqs.tre -t tips_to_keep.txt -o pruned.tre

**Prune a tree to remove the tips in tips_to_remove.txt. Note that the -n/--negate option must be passed for this functionality:**

::

	filter_tree.py -i rep_seqs.tre -t tips_to_keep.txt -o negated.tre -n

**Prune a tree to include only the tips found in the fasta file provided:**

::

	filter_tree.py -i rep_seqs.tre -f fast_f.fna -o pruned_fast.tre


