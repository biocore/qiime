.. _parallel_assign_taxonomy_uclust:

.. index:: parallel_assign_taxonomy_uclust.py

*parallel_assign_taxonomy_uclust.py* -- Parallel taxonomy assignment using the uclust consensus taxonomy assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `assign_taxonomy.py <./assign_taxonomy.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_assign_taxonomy_uclust.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to fasta file containing query sequences [REQUIRED]
	-o, `-`-output_dir
		Path to store output files [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-id_to_taxonomy_fp
		Full path to id_to_taxonomy mapping file [default: /Users/caporaso/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt]
	-r, `-`-reference_seqs_fp
		Ref seqs to search against. [default: /Users/caporaso/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta]
	`-`-min_consensus_fraction
		Minimum fraction of database hits that must have a specific taxonomic assignment to assign that taxonomy to a query [default: 0.51]
	`-`-similarity
		Minimum percent similarity to consider a database match a hit, expressed as a fraction between 0 and 1 [default: 0.9]
	`-`-uclust_max_accepts
		Number of database hits to consider when making an assignment [default: 3]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 1]
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

Mapping of sequence identifiers to taxonomy and quality information.


**Example:**

Assign taxonomy to all sequences in the input file (-i) using the uclust consensus taxonomy assigner and write the results (-o) to $PWD/uclust_assigned_taxonomy/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_assign_taxonomy_uclust.py -i $PWD/inseqs.fasta -o $PWD/uclust_assigned_taxonomy/


