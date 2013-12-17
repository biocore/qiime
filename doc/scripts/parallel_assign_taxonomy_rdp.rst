.. _parallel_assign_taxonomy_rdp:

.. index:: parallel_assign_taxonomy_rdp.py

*parallel_assign_taxonomy_rdp.py* -- Parallel taxonomy assignment using RDP
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `assign_taxonomy.py <./assign_taxonomy.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_assign_taxonomy_rdp.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to input_fasta_fp [REQUIRED]
	-o, `-`-output_dir
		Path to store output files [REQUIRED]
	
	**[OPTIONAL]**
		
	`-`-rdp_classifier_fp
		Full path to rdp classifier jar file [default: /Applications/rdp_classifier_2.2/rdp_classifier-2.2.jar]
	-c, `-`-confidence
		Minimum confidence to record an assignment [default: 0.8]
	-t, `-`-id_to_taxonomy_fp
		Full path to id_to_taxonomy mapping file [default: /Users/caporaso/data/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt]
	-r, `-`-reference_seqs_fp
		Ref seqs to rdp against. [default: /Users/caporaso/data/gg_13_5_otus/rep_set/97_otus.fasta]
	`-`-rdp_max_memory
		Maximum memory allocation, in MB, for Java virtual machine when using the rdp method.  Increase for large training sets [default: 4000]
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

Mapping of sequence identifiers to taxonomy and quality scores.


**Example:**

Assign taxonomy to all sequences in the input file (-i) using the RDP classifier and write the results (-o) to $PWD/rdp_assigned_taxonomy/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_assign_taxonomy_rdp.py -i $PWD/inseqs.fasta -o $PWD/rdp_assigned_taxonomy/


