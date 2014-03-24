.. _parallel_beta_diversity:

.. index:: parallel_beta_diversity.py

*parallel_beta_diversity.py* -- Parallel beta diversity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `beta_diversity.py <./beta_diversity.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_beta_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path, must be directory [REQUIRED]
	-o, `-`-output_path
		Output path, must be directory [REQUIRED]
	-m, `-`-metrics
		Metrics to use [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-tree_path
		Path to newick tree file, required for phylogenetic metrics [default: None]
	-N, `-`-beta_diversity_fp
		Full path to scripts/`beta_diversity.py <./beta_diversity.html>`_ [default: /python_software/Qiime/scripts/`beta_diversity.py <./beta_diversity.html>`_]
	-P, `-`-poller_fp
		Full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /python_software/Qiime/scripts/`poller.py <./poller.html>`_]
	-R, `-`-retain_temp_files
		Retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		Path to `cluster_jobs.py <./cluster_jobs.html>`_ script  [default: /python_software/Qiime/scripts/`start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		Full path to python executable [default: /opt/local/bin/python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 60]


**Output:**

The output of %prog is a folder containing text files, each a distance matrix between samples.


**Example:**

Apply the dist_unweighted_unifrac and the dist_weighted_unifrac metrics (-m) to all otu tables in /home/qiime_user/rare/ (-i) and write the resulting output files to /home/qiime_user/out/ (-o, will be created if it doesn't exist). Use the tree file /home/qiime_user/rep_set.tre (-t) when necessary.

::

	parallel_beta_diversity.py -i /home/qiime_user/rare/ -o /home/qiime_user/out -m dist_unweighted_unifrac,dist_weighted_unifrac -t /home/qiime_user/rep_set.tre


