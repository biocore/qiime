.. _core_qiime_analyses:

.. index:: core_qiime_analyses.py

*core_qiime_analyses.py* -- A workflow for running a core set of QIIME analyses.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script plugs several QIIME steps together to form a basic full data analysis workflow. The steps include quality filtering and demultiplexing sequences (optional), running the `pick_otus_through_otu_table.py <./pick_otus_through_otu_table.html>`_ workflow (pick otus and representative sequences, assign taxonomy, align representative sequences, build a tree, and build the OTU table), generating 2d and 3d beta diversity PCoA plots, generating alpha rarefaction plots, identifying OTUs that are differentially represented in different categories, and several additional analysis. Beta diversity calculations will be run both with and without an even sampling step, where the depth of sampling can either be passed to the script or QIIME will try to make a reasonable guess.


**Usage:** :file:`core_qiime_analyses.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fnas
		The input fasta file(s) -- comma-separated if more than one [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	-m, `-`-mapping_fp
		The mapping filepath [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters. [if omitted, default values will be used]
	-q, `-`-input_quals
		The 454 qual files. Comma-separated if more than one, and must correspond to the  order of the fasta files. Not relevant if passing  --suppress_split_libraries. [default: None]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-a, `-`-parallel
		Run in parallel where available. Specify number of jobs to start with -O or in the parameters file. [default: False]
	-e, `-`-seqs_per_sample
		Depth of coverage for diversity analyses that incorporate subsampling the OTU table to an equal number of sequences per sample. [default: determined automatically - bad choices can be made in some circumstances]
	`-`-even_sampling_keeps_all_samples
		If -e/--seqs_per_sample is not provided, chose the even sampling depth to force retaining all samples (rather then default which will choose a sampling depth which may favor keeping  more sequences by excluding some samples) [default: False]
	-t, `-`-reference_tree_fp
		Path to the tree file if one should be used. Relevant for closed-reference-based OTU picking methods only (i.e., uclust_ref -C and BLAST) [default: de novo tree will be used]
	-c, `-`-categories
		The metadata category or categories to compare (i.e., column headers in the mapping file) for the otu_category_significance, `supervised_learning.py <./supervised_learning.html>`_, and `cluster_quality.py <./cluster_quality.html>`_ steps. Pass a comma-separated list if more than one category [default: None; skip these steps]
	`-`-suppress_split_libraries
		Skip demultiplexing/quality filtering (i.e. split_libraries). This assumes that sequence identifiers are in post-split_libraries format (i.e., sampleID_seqID) [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 4]


**Output:**




Run serial analysis using a custom parameters file (-p), and guess the even sampling depth (no -e provided). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	core_qiime_analyses.py -i $PWD/Fasting_Example.fna -q $PWD/Fasting_Example.qual -o $PWD/FastingStudy_w_sl -m $PWD/Fasting_Map.txt -c BarcodeSequence -p $PWD/params.txt

Run serial analysis using a custom parameters file (-p), and guess the even sampling depth (no -e provided). Skip split libraries by starting with already demultiplexed sequences. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	core_qiime_analyses.py -i $PWD/seqs.fna -o $PWD/FastingStudy -m $PWD/Fasting_Map.txt -c BarcodeSequence --suppress_split_libraries -p $PWD/params.txt


