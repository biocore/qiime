.. _parallel_blast:

.. index:: parallel_blast

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
	-r, `-`-refseqs_path
		Path to fasta sequences to search against [REQUIRED]
	-o, `-`-output_dir
		name of output directory for blast jobs [REQUIRED]
	
	**[OPTIONAL]**
		
	-e, `-`-e_value
		E-value threshold for blasts [default: 1e-30]
	-n, `-`-num_hits
		number of hits per query for blast results [default: 1]
	-w, `-`-word_size
		word size for blast searches [default: 30]
	-D, `-`-suppress_format_blastdb
		supress format of blastdb [default: False]
	-a, `-`-blastmat_dir
		full path to directory containing blastmat file [default: /Users/Jesse/blast-2.2.21/data]
	-b, `-`-blastall_fp
		Path to blastall [default: /Users/Jesse/blast-2.2.21/bin/blastall]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 24]
	-P, `-`-poller_fp
		full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /Users/Jesse/Qiime/qiime/parallel/`poller.py <./poller.html>`_]
	-R, `-`-retain_temp_files
		retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		path to `cluster_jobs.py <./cluster_jobs.html>`_ script  [default: /software/scripts/`cluster_jobs.py <./cluster_jobs.html>`_]
	-W, `-`-suppress_polling
		suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		full path to python executable [default: /usr/local/bin/python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 60]


**Output:**

 


**Example**

Split 10_seq.fasta (-i) into three fasta files (-O) and blast each against blast database created from 1000_seq.fasta (-d)

::

	parallel_blast.py -i 10_seq.fasta -d 1000_seq.fasta -O 3 -o bla_out


