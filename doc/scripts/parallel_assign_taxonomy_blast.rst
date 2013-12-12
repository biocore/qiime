.. _parallel_assign_taxonomy_blast:

.. index:: parallel_assign_taxonomy_blast.py

*parallel_assign_taxonomy_blast.py* -- Parallel taxonomy assignment using BLAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `assign_taxonomy.py <./assign_taxonomy.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_assign_taxonomy_blast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to input_fasta_fp [REQUIRED]
	-o, `-`-output_dir
		Full path to store output files [REQUIRED]
	
	**[OPTIONAL]**
		
	-r, `-`-reference_seqs_fp
		Ref seqs to blast against.  Must provide either --blast_db or --reference_seqs_db for assignment with blast [default: /Users/caporaso/data/gg_13_5_otus/rep_set/97_otus.fasta]
	-b, `-`-blast_db
		Database to blast against.  Must provide either --blast_db or --reference_seqs_db for assignment with blast [default: None]
	-e, `-`-e_value
		Maximum e-value to record an assignment, only used for blast method [default: 0.001]
	-B, `-`-blastmat_dir
		Full path to directory containing blastmat file [default: /Applications/blast-2.2.22/data/]
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
	-t, `-`-id_to_taxonomy_fp
		Full path to id_to_taxonomy mapping file [default: /Users/caporaso/data/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt]


**Output:**

Mapping of sequence identifiers to taxonomy and quality scores.


**Example:**

Assign taxonomy to all sequences in the input file (-i) using BLAST with the id to taxonomy mapping file (-t) and reference sequences file (-r), and write the results (-o) to $PWD/blast_assigned_taxonomy/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_assign_taxonomy_blast.py -i $PWD/inseqs.fasta -t $PWD/id_to_tax.txt -r $PWD/refseqs.fasta -o $PWD/blast_assigned_taxonomy/


