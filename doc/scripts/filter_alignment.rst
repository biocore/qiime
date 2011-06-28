.. _filter_alignment:

.. index:: filter_alignment.py

*filter_alignment.py* -- Filter sequence alignment by removing highly variable regions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script should be applied to generate a useful tree when aligning against a template alignment (e.g., with PyNAST). This script will remove positions which are gaps in every sequence (common for PyNAST, as typical sequences cover only 200-400 bases, and they are being aligned against the full 16S gene). Additionally, the user can supply a lanemask file, that defines which positions should included when building the tree, and which should be ignored. Typically, this will differentiate between non-conserved positions, which are uninformative for tree building, and conserved positions which are informative for tree building. FILTERING ALIGNMENTS WHICH WERE BUILD WITH PYNAST AGAINST THE GREENGENES CORE SET ALIGNMENT SHOULD BE CONSIDERED AN ESSENTIAL STEP.


**Usage:** :file:`filter_alignment.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_file
		The input directory 
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		The output directory [default: .]
	-m, `-`-lane_mask_fp
		Path to lanemask file [default: None]
	-s, `-`-suppress_lane_mask_filter
		Suppress lane mask filtering (necessary to turn off lane-mask-based filtering when a qiime_config default is  provided for --lane_mask_fp) [default: False]
	-g, `-`-allowed_gap_frac
		Gap filter threshold, filters positions which are gaps in > allowed_gap_frac of the sequences [default: 0.999999]
	-r, `-`-remove_outliers
		Remove seqs very dissimilar to the alignment consensus (see --threshold).  [default: False]
	-t, `-`-threshold
		With -r, remove seqs whose dissimilarity to the consensus sequence is approximately > x standard devaitions above the mean of the sequences [default: 3.0]
	-e, `-`-entropy_threshold
		Sets percent threshold for removing base positions with the highest entropy.  For example, if 0.10 were specified, the top 10% most entropic base positions would be filtered.  If this value is used, any lane mask supplied will be ignored.  Entropy filtered occurs after gap filtering.    [default: None]


**Output:**

The output of `filter_alignment.py <./filter_alignment.html>`_ consists of a single FASTA file, which ends with "pfiltered.fasta", where the "p" stands for positional filtering of the columns.


**Examples:**

As a simple example of this script, the user can use the following command, which consists of an input FASTA file (i.e. resulting file from `align_seqs.py <./align_seqs.html>`_), lanemask template file and the output directory "filtered_alignment/":

::

	filter_alignment.py -i repr_set_seqs_aligned.fna -m lanemask_template -o filtered_alignment/

Alternatively, if the user would like to use a different gap fraction threshold ("-g"), they can use the following command:

::

	filter_alignment.py -i repr_set_seqs_aligned.fna -m lanemask_template -o filtered_alignment/ -g 0.95


