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
		full path to input_fasta_fp [REQUIRED]
	-o, `-`-output_dir
		path to store output files [REQUIRED]
	
	**[OPTIONAL]**
		
	-r, `-`-rdp_classifier_fp
		full path to rdp classifier jar file [default: /Users/Jesse/rdp_classifier/rdp_classifier-2.0.jar]
	-c, `-`-confidence
		Minimum confidence to record an assignment [default: 0.8]
	-N, `-`-assign_taxonomy_fp
		full path to scripts/`assign_taxonomy.py <./assign_taxonomy.html>`_ [default: /Users/Jesse/Qiime/scripts/`assign_taxonomy.py <./assign_taxonomy.html>`_]
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

Mapping of sequence identifiers to taxonomy and quality scores.


**Example**

Split the input file (-i) into five jobs (-O) start them,and write the results (-o) to out/. BE SURE TO SPECIFY FULL PATHS!

::

	parallel_assign_taxonomy_rdp -O 5 -i inseqs.fasta -o out/


