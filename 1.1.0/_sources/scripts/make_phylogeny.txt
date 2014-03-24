.. _make_phylogeny:

.. index:: make_phylogeny.py

*make_phylogeny.py* -- Make Phylogeny
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Many downstream analyses require that the phylogenetic tree relating the OTUs in a study be present. The script `make_phylogeny.py <./make_phylogeny.html>`_ produces this tree from a multiple sequence alignment. Trees are constructed with a set of sequences representative of the OTUs, by default using FastTree (Price, Dehal, & Arkin, 2009).


**Usage:** :file:`make_phylogeny.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		Path to read input alignment
	
	**[OPTIONAL]**
		
	-t, `-`-tree_method
		Method for tree building. Valid choices are: clearcut, clustalw, raxml, fasttree_v1, fasttree, muscle [default: fasttree]
	-o, `-`-result_fp
		Path to store result file [default: <input_sequences_filename>.tre]
	-l, `-`-log_fp
		Path to store log file [default: No log file created.]
	-r, `-`-root_method
		Method for choosing root of phylo tree  Valid choices are: midpoint, tree_method_default [default: tree_method_default]


**Output:**

The result of `make_phylogeny.py <./make_phylogeny.html>`_ consists of a newick formatted tree file (.tre) and optionally a log file. The tree file is formatted using the Newick format and this file can be viewed using most tree visualization tools, such as TopiaryTool, FigTree, etc.


**Examples:**

A simple example of `make_phylogeny.py <./make_phylogeny.html>`_ is shown by the following command, where we use the default tree building method (fasttree) and write the file to the current working directory without a log file:

::

	make_phylogeny.py -i repr_set_seqs_aligned_pfiltered.fasta -o rep_phylo.tre

Alternatively, if the user would prefer using another tree building method (i.e. clearcut (Sheneman, Evans, & Foster, 2006)), then they could use the following command:

::

	make_phylogeny.py -i repr_set_seqs_aligned_pfiltered.fasta -t clearcut

Note: For whichever method used, the 3rd party program must be properly installed on the user's computer.


