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
		Path to read input fasta alignment, only first word in defline will be considered
	
	**[OPTIONAL]**
		
	-t, `-`-tree_method
		Method for tree building. Valid choices are: clearcut, clustalw, fasttree_v1, fasttree, raxml_v730, muscle [default: fasttree]
	-o, `-`-result_fp
		Path to store result file [default: <input_sequences_filename>.tre]
	-l, `-`-log_fp
		Path to store log file [default: No log file created.]
	-r, `-`-root_method
		Method for choosing root of phylo tree  Valid choices are: midpoint, tree_method_default [default: tree_method_default]


**Output:**

The result of `make_phylogeny.py <./make_phylogeny.html>`_ consists of a newick formatted tree file (.tre) and optionally a log file. The tree file is formatted using the Newick format and this file can be viewed using most tree visualization tools, such as TopiaryTool, FigTree, etc.

The tips of the tree are the first word from the input sequences from the fasta file, e.g.: '>101 PC.481_71 RC:1..220' is represented in the tree as '101'.


**Examples:**

A simple example of `make_phylogeny.py <./make_phylogeny.html>`_ is shown by the following command, where we use the default tree building method (fasttree) and write the file to the current working directory without a log file:

::

	make_phylogeny.py -i $PWD/aligned.fasta -o $PWD/rep_phylo.tre

Alternatively, if the user would prefer using another tree building method (i.e. clearcut (Sheneman, Evans, & Foster, 2006)), then they could use the following command:

::

	make_phylogeny.py -i $PWD/aligned.fasta -t clearcut


