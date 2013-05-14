.. _parallel_blast:

.. index:: parallel_blast.py

*parallel_blast.py* -- Parallel BLAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script for performing blast while making use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_blast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-infile_path
		Path of sequences to use as queries [REQUIRED]
	-o, `-`-output_dir
		Name of output directory for blast jobs [REQUIRED]
	
	**[OPTIONAL]**
		
	-c, `-`-disable_low_complexity_filter
		Disable filtering of low-complexity sequences (i.e., -F F is passed to blast) [default: False]
	-e, `-`-e_value
		E-value threshold for blasts [default: 1e-30]
	-n, `-`-num_hits
		Number of hits per query for blast results [default: 1]
	-w, `-`-word_size
		Word size for blast searches [default: 30]
	-a, `-`-blastmat_dir
		Full path to directory containing blastmat file [default: /Applications/blast-2.2.22/data/]
	-r, `-`-refseqs_path
		Path to fasta sequences to search against. Required if -b is not provided.
	-b, `-`-blast_db
		Name of pre-formatted BLAST database. Required if -r is not provided.
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

BLAST $PWD/inseqs.fasta (-i) against a blast database created from $PWD/refseqs.fasta (-r). Store the results in $PWD/blast_out/ (-o). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_blast.py -i $PWD/inseqs.fasta -r $PWD/refseqs.fasta -o $PWD/blast_out/ -e 0.001


