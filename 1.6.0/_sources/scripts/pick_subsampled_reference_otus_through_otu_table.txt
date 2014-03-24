.. _pick_subsampled_reference_otus_through_otu_table:

.. index:: pick_subsampled_reference_otus_through_otu_table.py

*pick_subsampled_reference_otus_through_otu_table.py* -- 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`pick_subsampled_reference_otus_through_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		The input sequences filepath or comma-separated list of filepaths
	-r, `-`-reference_fp
		The reference sequences
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	`-`-prefilter_refseqs_fp
		The reference sequences to use for the prefilter, if different from the reference sequecnces to use for the OTU picking [default: same as passed for --reference_fp]
	-n, `-`-new_ref_set_id
		Unique identifier for OTUs that get created in this ref set (this is useful to support combining of reference sets) [default:New]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 4]
	-s, `-`-percent_subsample
		Percent of failure sequences to include in the subsample to cluster de novo (larger numbers should give more comprehensive ,results but will be slower) [default:0.001]
	`-`-prefilter_percent_id
		Sequences are pre-clustered at this percent id against the reference and any reads which fail to hit are discarded (a quality filter); pass 0.0 to disable [default:0.60]
	`-`-step1_otu_map_fp
		Reference OTU picking OTU map  (to avoid rebuilding if one has already been built)
	`-`-step1_failures_fasta_fp
		Reference OTU picking failures fasta filepath  (to avoid rebuilding if one has already been built)
	`-`-suppress_step4
		Suppress the final de novo OTU picking step  (may be necessary for extremely large data sets) [default: False]
	`-`-min_otu_size
		The minimum otu size (in number of sequences) to retain the otu [default: 2]
	`-`-suppress_taxonomy_assignment
		Skip the taxonomy assignment step, resulting in an OTU table without taxonomy [default: False]
	`-`-suppress_align_and_tree
		Skip the sequence alignment and tree-building steps [default: False]


**Output:**




Run the subsampled open-reference OTU picking workflow on seqs1.fna using refseqs.fna as the reference collection. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_subsampled_reference_otus_through_otu_table.py -i $PWD/seqs1.fna -r $PWD/refseqs.fna -o $PWD/ucrss/ -s 0.1 -p $PWD/ucrss_params.txt

Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_subsampled_reference_otus_through_otu_table.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter/ -s 0.1 -p $PWD/ucrss_params.txt

Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection. This is useful if you're working with marker genes that do not result in useful alignment (e.g., fungal ITS). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_subsampled_reference_otus_through_otu_table.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter_no_tree/ -s 0.1 -p $PWD/ucrss_params.txt --suppress_align_and_tree

Run the subsampled open-reference OTU picking workflow in iterative mode on seqs1.fna and seqs2.fna using refseqs.fna as the initial reference collection, suppressing assignment of taxonomy. This is useful if you're working with a reference collection without associated taxonomy. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/

::

	pick_subsampled_reference_otus_through_otu_table.py -i $PWD/seqs1.fna,$PWD/seqs2.fna -r $PWD/refseqs.fna -o $PWD/ucrss_iter_no_tax/ -s 0.1 -p $PWD/ucrss_params.txt --suppress_taxonomy_assignment


