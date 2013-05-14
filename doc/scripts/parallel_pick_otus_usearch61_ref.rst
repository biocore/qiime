.. _parallel_pick_otus_usearch61_ref:

.. index:: parallel_pick_otus_usearch61_ref.py

*parallel_pick_otus_usearch61_ref.py* -- Parallel pick otus using usearch_ref
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script works like the `pick_otus.py <./pick_otus.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_pick_otus_usearch61_ref.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to input_fasta_fp
	-o, `-`-output_dir
		Path to store output files
	-r, `-`-refseqs_fp
		Full path to reference collection
	
	**[OPTIONAL]**
		
	-s, `-`-similarity
		Sequence similarity threshold [default: 0.97]
	-z, `-`-enable_rev_strand_match
		Enable reverse strand matching for uclust, uclust_ref, usearch, usearch_ref, usearch61, or usearch61_ref otu picking, will double the amount of memory used. [default: False]
	`-`-max_accepts
		Max_accepts value to uclust, uclust_ref, usearch61, and usearch61_ref.  By default, will use value suggested by method (uclust: 20, usearch61: 1) [default: default]
	`-`-max_rejects
		Max_rejects value for uclust, uclust_ref, usearch61, and usearch61_ref.  With default settings, will use value recommended by clustering method used (uclust: 500, usearch61: 8 for usearch_fast_cluster option, 32 for reference and smallmem options) [default: default]
	`-`-word_length
		Word length value for uclust, uclust_ref, and usearch, usearch_ref, usearch61, and usearch61_ref. With default setting, will use the setting recommended by the method (uclust: 12, usearch: 64, usearch61: 8).  int value can be supplied to override this setting. [default: default]
	`-`-minlen
		Minimum length of sequence allowed for usearch, usearch_ref, usearch61, and usearch61_ref. [default: 64]
	`-`-usearch_fast_cluster
		Use fast clustering option for usearch or usearch61_ref with new clusters.  --enable_rev_strand_match can not be enabled with this option, and the only valid option for usearch61_sort_method is 'length'.  This option uses more memory than the default option for de novo clustering. [default: False]
	`-`-usearch61_sort_method
		Sorting method for usearch61 and usearch61_ref.  Valid options are abundance, length, or None.  If the --usearch_fast_cluster option is enabled, the only sorting method allowed in length. [default: abundance]
	`-`-sizeorder
		Enable size based preference in clustering with usearch61. Requires that --usearch61_sort_method be abundance. [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]
	-R, `-`-retain_temp_files
		Retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		Path to cluster jobs script (defined in qiime_config)  [default: `start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 1]


**Output:**




**Example:**

Pick OTUs by searching $PWD/inseqs.fasta against $PWD/refseqs.fasta with reference-based usearch and write the output to the $PWD/usearch_ref_otus/ directory. This is a closed-reference OTU picking process. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_pick_otus_usearch61_ref.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/usearch_ref_otus/


