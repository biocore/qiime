.. _pick_otus_through_otu_table:

.. index:: pick_otus_through_otu_table.py

*pick_otus_through_otu_table.py* -- A workflow script for picking OTUs through building OTU tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a sequence file and performs all processing steps through building the OTU table.

REQUIRED: You must add values for the following parameters in a custom parameters file:
 align_seqs:template_fp
 filter_alignment:lane_mask_fp 
 
These are the values that you would typically pass as --template_fp to `align_seqs.py <./align_seqs.html>`_ and lane_mask_fp to `filter_alignment.py <./filter_alignment.html>`_, respectively.





**Usage:** :file:`pick_otus_through_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input fasta file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	-p, `-`-parameter_fp
		Path to the parameter file [REQUIRED]
	
	**[OPTIONAL]**
		
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-m, `-`-mapping_fp
		The mapping filepath [REQUIRED for denoising]
	-s, `-`-sff_fp
		The sff file [REQUIRED for denoising]


**Output:**

This script will produce a set of cluster centroids (as a FASTA file) and a cluster mapping file (from `denoise.py <./denoise.html>`_ if sff.txt and mapping file were provided), an OTU mapping file (`pick_otus.py <./pick_otus.html>`_), a representative set of sequences (FASTA file from `pick_rep_set.py <./pick_rep_set.html>`_), a sequence alignment file (FASTA file from `align_seqs.py <./align_seqs.html>`_), taxonomy assignment file (from `assign_taxonomy.py <./assign_taxonomy.html>`_), a filtered sequence alignment (from `filter_alignment.py <./filter_alignment.html>`_), a phylogenetic tree (Newick file from `make_phylogeny.py <./make_phylogeny.html>`_) and an OTU table (from `make_otu_table.py <./make_otu_table.html>`_).


**Simple example:**

The following command will start an analysis on inseq1.fasta (-i), which is a post-split_libraries fasta file. The sequence identifiers in this file should be of the form <sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data preparation, are applied.

1. Pick OTUs with uclust at similarity of 0.97;

2. Pick a representative set with the most_abundant method;

3. Align the representative set with PyNAST (REQUIRED: SPECIFY TEMPLATE ALIGNMENT with align_seqs:template_fp in the parameters file);

4. Assign taxonomy with RDP classifier;

5. Filter the alignment prior to tree building - remove positions which are all gaps, and specified as 0 in the lanemask (REQUIRED: SPECIFY LANEMASK with filter_alignment:lane_mask_fp in the parameters file);

6. Build a phylogenetic tree with FastTree;

7. Build an OTU table.

All output files will be written to the directory specified by -o, and 
subdirectories as appropriate.


::

	pick_otus_through_otu_table.py -i inseqs1.fasta -o wf1/ -p custom_parameters.txt

**Simple example with denoising:**

This command will do the same steps as the previous example and additionally denoise the data set prior to OTU picking. Only flowgrams in the input sff.txt file that have a matching identifier in inseqs1.fasta are considered here, the rest is discarded.

1. Denoise flowgrams in inseqs1.sff.txt;

2. Pick OTUs with uclust at similarity of 0.97;

3. Pick a representative set with the most_abundant method;

4. Align the representative set with PyNAST (REQUIRED: SPECIFY TEMPLATE ALIGNMENT with align_seqs:template_fp in the parameters file);

5. Assign taxonomy with RDP classifier;

6. Filter the alignment prior to tree building - remove positions which are all gaps, and specified as 0 in the lanemask (REQUIRED: SPECIFY LANEMASK with filter_alignment:lane_mask_fp in the parameters file);

7. Build a phylogenetic tree with FastTree;

8. Build an OTU table.

All output files will be written to the directory specified by -o, and 
subdirectories as appropriate.


::

	pick_otus_through_otu_table.py -f inseqs1.sff.txt -m metadata_mapping.txt -i inseqs1.fasta -o wf2/ -p custom_parameters.txt


